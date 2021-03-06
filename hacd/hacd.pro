QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

include($$PWD/HACD_Lib/HACD_Lib.pri)
INCLUDEPATH += $$PWD/HACD_Lib/inc/

HEADER_FILES = $$files($$PWD/*.h, false)
SOURCE_FILES = $$files($$PWD/*.cpp, false)

HEADERS += $${HEADER_FILES}
SOURCES += $${SOURCE_FILES}

DEFINES += PRO_PATH=\"\\\"$$PWD\\\"\"


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
