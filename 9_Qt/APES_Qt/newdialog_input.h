#ifndef NEWDIALOG_H
#define NEWDIALOG_H

#include <QDialog>
#include <QTreeWidget>
#include <QJsonObject>

class NewDialog : public QDialog
{
    Q_OBJECT

public:
    explicit NewDialog(QWidget *parent = nullptr);
    void setJson(QString infilePath);

private:
    QTreeWidget *treeWidget;
    QString filePath;
    QJsonObject currentJsonObject;
private slots:
    void onItemChanged(QTreeWidgetItem *item, int column);
    void updateJson(const QString &path, const QString &newValue,QJsonObject &currentJsonObject);
    void onUpdateButtonClicked();
signals:
    void confirmClicked();
};

#endif // NEWDIALOG_H
