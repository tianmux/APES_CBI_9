#include "newdialog_input.h"
#include <QVBoxLayout>
#include <QJsonValue>
#include <QJsonObject>
#include <QJsonArray>
#include <QFile>
#include <QJsonDocument>
#include <QString>
#include <QPushButton>

QJsonObject readJsonFile(const QString &fileName) {
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly)) {
        qWarning("Failed to open file");
        return QJsonObject();
    }

    QByteArray data = file.readAll();
    QJsonDocument doc(QJsonDocument::fromJson(data));
    return doc.object();
}

bool writeJsonToFile(const QJsonObject &jsonObject, const QString &filePath) {
    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly)) {
        qWarning("Couldn't open file for writing.");
        return false;
    }

    QJsonDocument doc(jsonObject);
    file.write(doc.toJson(QJsonDocument::Indented)); // Use QJsonDocument::Compact for non-indented
    file.close();

    return true;
}
void addItems(QTreeWidgetItem *parentItem, const QJsonValue &value) {
    if (value.isObject()) {
        QJsonObject obj = value.toObject();
        for (auto it = obj.begin(); it != obj.end(); ++it) {
            QTreeWidgetItem *childItem = new QTreeWidgetItem();
            childItem->setText(0, it.key());
            addItems(childItem, it.value());
            childItem->setFlags(childItem->flags() | Qt::ItemIsEditable); // Make value editable
            parentItem->addChild(childItem);

        }
    } else if (value.isArray()) {
        QJsonArray array = value.toArray();
        for (int i = 0; i < array.size(); ++i) {
            QTreeWidgetItem *childItem = new QTreeWidgetItem();
            childItem->setText(0, QString("[%1]").arg(i));
            addItems(childItem, array[i]);
            childItem->setFlags(childItem->flags() | Qt::ItemIsEditable); // Make value editable
            parentItem->addChild(childItem);
        }
    } else {
        parentItem->setText(1, value.toVariant().toString());
    }
}

void populateTreeWidget(QTreeWidget *treeWidget, const QJsonObject &jsonObject) {
    treeWidget->clear();
    treeWidget->setColumnCount(2);
    treeWidget->setHeaderLabels(QStringList() << "Property" << "Value");

    for (auto it = jsonObject.begin(); it != jsonObject.end(); ++it) {
        QTreeWidgetItem *item = new QTreeWidgetItem(treeWidget);
        item->setText(0, it.key());
        addItems(item, it.value());
    }
}
NewDialog::NewDialog(QWidget *parent) : QDialog(parent), treeWidget(new QTreeWidget(this)) {
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->addWidget(treeWidget);
    QPushButton *updateButton = new QPushButton("Confirm", this);

    layout->addWidget(updateButton);
    updateButton->setGeometry(1000,100,100,25);

    connect(updateButton, &QPushButton::clicked, this, &NewDialog::onUpdateButtonClicked);
    connect(treeWidget, &QTreeWidget::itemChanged, this, &NewDialog::onItemChanged);

    treeWidget->setColumnCount(2);
    treeWidget->setHeaderLabels(QStringList() << "Property" << "Value");
    // Additional setup for treeWidget if needed
}

void NewDialog::onUpdateButtonClicked() {
    if (writeJsonToFile(currentJsonObject, this->filePath)) {
        emit confirmClicked();
    } else {
        qDebug() << "Failed to write JSON data to" << filePath;
    }
}

void NewDialog::setJson(QString infilePath) {
    this->filePath = infilePath;
    // Clear existing items
    treeWidget->clear();
    treeWidget->setGeometry(1000,100,800,1000);
    this->currentJsonObject=readJsonFile(filePath);
    // Function to populate treeWidget with JSON data
    // Similar to the populateTreeWidget function shown earlier
    populateTreeWidget(treeWidget, currentJsonObject);
    treeWidget->expandAll();
    int numberOfColumns = treeWidget->columnCount();
    int treeWidgetWidth = treeWidget->width();

    for (int i = 0; i < numberOfColumns; ++i) {
        treeWidget->setColumnWidth(i, treeWidgetWidth / numberOfColumns);
    }
}

void NewDialog::onItemChanged(QTreeWidgetItem *item, int column) {
    if (column != 1) {
        // Only values in the second column are editable
        return;
    }

    // Build the path to the JSON value that this tree item represents
    QStringList pathParts;
    QTreeWidgetItem *currentItem = item;
    while (currentItem) {
        pathParts.prepend(currentItem->text(0)); // Add the key to the path
        currentItem = currentItem->parent();    // Move up to the parent item
    }

    // Combine path parts into a single string
    QString path = pathParts.join(".");

    // Update the JSON object with the new value
    updateJson(path, item->text(1),this->currentJsonObject);
}

void NewDialog::updateJson(const QString &path, const QString &newValue,QJsonObject &currentJsonObject) {
    QStringList keys = path.split('.', Qt::SkipEmptyParts);
    QJsonValueRef ref = currentJsonObject[keys[0]];

    // Traverse the JSON using the path, but stop at the parent of the target element
    for (int i = 1; i < keys.size() - 1; ++i) {
        if (ref.isObject()) {
            ref = ref.toObject().value(keys[i]);
        } else if (ref.isArray()) {
            int index = keys[i].toInt();
            ref = ref.toArray().at(index);
        }
    }

    // Convert newValue to the appropriate type if necessary
    QJsonValue updatedValue;
    if (ref.isObject()) {
        QJsonObject obj = ref.toObject();
        QJsonValue oldValue = obj.value(keys.last());

        if (oldValue.isDouble()) {
            bool ok;
            double doubleValue = newValue.toDouble(&ok);
            updatedValue = ok ? QJsonValue(doubleValue) : QJsonValue(newValue);
        } else if (oldValue.isBool()) {
            updatedValue = (newValue.toLower() == "true");
        } else {
            updatedValue = QJsonValue(newValue);
        }

        obj[keys.last()] = updatedValue;
        ref = obj; // Important: update the parent reference
    } else if (ref.isArray()) {
        QJsonArray arr = ref.toArray();
        int index = keys.last().toInt(); // Assuming the last key is an array index
        updatedValue = QJsonValue(newValue); // Simplified: treating as string

        arr[index] = updatedValue;
        ref = arr; // Important: update the parent reference
    }
}
