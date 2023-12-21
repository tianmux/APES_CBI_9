#include "chartwindow.h"

ChartWindow::ChartWindow(QWidget *parent) :
    QWidget(parent),
    chartView(new QChartView(this))
{
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->addWidget(chartView);
}

void ChartWindow::setupChart(QAbstractSeries *series, const QString &title)
{
    QChart *chart = new QChart();
    chart->addSeries(series);
    chart->setTitle(title);
    chart->createDefaultAxes();
    chartView->setChart(chart);
}
