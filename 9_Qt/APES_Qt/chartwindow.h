#ifndef CHARTWINDOW_H
#define CHARTWINDOW_H

#include <QWidget>
#include <QtCharts>

class ChartWindow : public QWidget
{
    Q_OBJECT

public:
    explicit ChartWindow(QWidget *parent = nullptr);
    void setupChart(QAbstractSeries *series, const QString &title);

private:
    QChartView *chartView;
};

#endif // CHARTWINDOW_H
