#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include <QLineEdit>
#include "ui_mainwindow.h"  // Include the generated UI header
#include "InputData.h"
#include "widgetstream.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;

}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    MainWindow(QWidget *parent = nullptr);
    InputData inputData;
    WidgetStream *myStream;
    ~MainWindow();

private slots:
    int on_pushButton_clicked();

    int on_startSim_clicked();

    void on_plotFFT_clicked();
    void enableStartSimButton();
private:
    Ui::MainWindow *ui;
    QString filePath;
};
#endif // MAINWINDOW_H
