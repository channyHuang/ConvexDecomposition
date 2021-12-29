#include "widget.h"

#include <QApplication>

#include <Qdir>
#include <QDebug>


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Widget w;
    w.show();
    return a.exec();
}
