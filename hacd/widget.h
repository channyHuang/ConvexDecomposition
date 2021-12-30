#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include <QFileDialog>
#include <QLineEdit>
#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QComboBox>

class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent = nullptr);
    ~Widget();

    void initConnection();

public slots:
    void selectPath();
    void action();

private:
    QLineEdit *pathLabel;
    QPushButton *pathBtn, *singleDirBtn;
    QComboBox *comboBox;
};
#endif // WIDGET_H
