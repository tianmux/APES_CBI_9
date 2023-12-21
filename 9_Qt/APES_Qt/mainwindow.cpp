#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFileDialog>
#include <QPushButton>
#include <QLineEdit>
#include <iostream>
#include "Beam.h"
#include "Cavity.h"
#include "Ring.h"
#include "PhysicsModule.h"
#include "InputData.h"
#include <fftw3.h>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include "chartwindow.h"
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QJsonValue>
#include <QJsonObject>
#include <QJsonArray>
#include <QFile>
#include <QJsonDocument>
#include "newdialog_input.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->myStream = new WidgetStream(std::cout,ui->textEdit_output);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::enableStartSimButton() {
    ui->startSim->setEnabled(true); // Assuming the button is named startSimButton
}
int  MainWindow::on_pushButton_clicked()
{
    this->filePath = QFileDialog::getOpenFileName(this, "Select File");
    if (!filePath.isEmpty()) {
        ui->lineEdit->setText(filePath);
        NewDialog *dialog = new NewDialog(this);
        connect(dialog, &NewDialog::confirmClicked, this, &MainWindow::enableStartSimButton);
        dialog->setJson(filePath); // jsonObject is your QJsonObject
        dialog->resize(800,600);
        dialog->setGeometry(1000,100,800,1000);
        dialog->show(); // Show the dialog as modal
        //populateTreeWidget(ui->treeWidget, jsonObject);
    }
    return 0;
}

int MainWindow::on_startSim_clicked()
{

    if (!filePath.isEmpty()) {
        ui->lineEdit->setText(filePath);
        // Load data from a test JSON file
        std::string inputFilePath = filePath.toStdString();
        try {
            inputData.loadDataFromFile(inputFilePath);
        } catch (const std::exception& e) {
            std::cerr << "Failed to load input data: " << e.what() << std::endl;
            ui->readInput_status->setText("Read input file fail.");
            ui->startSim->setEnabled(false);
            return 1;
        }
        ui->readInput_status->setText("Read input file successfully.");
        ui->startSim->setEnabled(true);
        //QJsonObject jsonObject = readJsonFile(filePath);
    }
    // Initialize the Beam object
    Beam beam(inputData);

    // Initialize the Cavity object
    Cavity cavity(inputData,beam);

    // Initialize the Ring object
    Ring ring(inputData);

    // Print the loaded data
    std::cout<< "Printing the loaded data" << std::endl;
    inputData.printInputData();
    PhysicsModule physicsModule;
    // Perform one turn tracking
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < inputData.latticeProps.nTrack; i++){
        if(i%int(inputData.latticeProps.nTrack/10)==0){
            end = std::chrono::high_resolution_clock::now();
            std::cout<<"Turn: "<< i << std::endl;
            std::cout << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;
        }
        physicsModule.oneTurn(inputData, cavity, ring, beam);
        //beam.storeParCentroids(i);
        //cavity.printCavityProps();
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken to finish full simulation: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " seconds" << std::endl;
    // write the centroids to file
    beam.writeCentroidsToFile(inputData.generalProps.projectName+"_centroid.txt",0,1);

    return 0;
}


void MainWindow::on_plotFFT_clicked()
{
    QVector<double> timeData;
    QVector<double> gammaData;
    QString filePath = QString::fromStdString(inputData.generalProps.projectName+"_centroid.txt");
    QFile file(filePath);
    if(file.open(QIODevice::ReadOnly)) {
        QTextStream in(&file);
        std::cout<<"Reading data..."<<std::endl;
        while(!in.atEnd()) {
            QString line = in.readLine();
            QStringList fields = line.split(' ');
            if(fields.size() >= 2) {
                timeData.append(fields[0].toDouble());
                gammaData.append(fields[1].toDouble());
            }
        }
        file.close();
        std::cout<<"Read the data."<<std::endl;
        std::cout<<timeData[0]<<std::endl;
        // Read the data
        fftw_complex *inData, *outData;
        fftw_plan p;
        QString text = ui->FFT_P->text();

        int N = text.toInt();//timeData.length();
        inData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N); // N is the number of points
        outData = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

        // Populate 'inData' with your data
        for(int i = 0; i < N; i++) {
            inData[i][0] = gammaData[i]; // Real part
            inData[i][1] = 0.0; // Imaginary part
        }

        std::cout<<"Put the data in raw vector."<<std::endl;
        std::cout<<"The length of the data is "<<N<<std::endl;
        p = fftw_plan_dft_1d(N, inData, outData, FFTW_FORWARD, FFTW_ESTIMATE);

        fftw_execute(p); // Execute the FFT

        std::cout<<"Done FFT."<<std::endl;

        // Plot
        QVector<QPointF> fftPoints; // For storing processed FFT data

        // Assuming 'out' is your FFTW output and 'N' is the number of points
        double maxMagnitude = 0.0;
        int maxIndex = 0;

        for (int i = 0; i < N; ++i) {
            double realPart = outData[i][0];
            double imagPart = outData[i][1];
            double magnitude = sqrt(realPart * realPart + imagPart * imagPart);

            fftPoints.append(QPointF(i, magnitude)); // For example, plotting magnitude against frequency bin 'i'

            if (magnitude > maxMagnitude) {
                maxMagnitude = magnitude;
                maxIndex = i;
            }
        }
        QLineSeries *series = new QLineSeries();
        for (const QPointF &point : fftPoints) {
            series->append(point);
        }
        ChartWindow *chartWindow = new ChartWindow();
        chartWindow->resize(800,600);
        chartWindow->setupChart(series, "FFT Results");
        chartWindow->setGeometry(100,600,800,400);
        chartWindow->show();

        // Find the peaks
        std::vector<int> peakIndices;

        double threshold = maxMagnitude*0.9; // Set a threshold value

        for (int i = 1; i < N - 1; ++i) {
            if (fftPoints[i].y() > threshold && fftPoints[i].y() > fftPoints[i-1].y()&& fftPoints[i].y() > fftPoints[i+1].y()) {
                peakIndices.push_back(i); // Store the index of the peak
            }
        }
        for (int idx : peakIndices) {
            double frequency = idx*inputData.fRev/N ;
            std::cout << "Peak at frequency: " << frequency << " Hz" << std::endl;
        }
    }
    else{
        std::cout<<"Could not find the ouput file."<<std::endl;
    }
}

