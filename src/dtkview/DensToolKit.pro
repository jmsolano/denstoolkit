#-------------------------------------------------
#
# Project created by QtCreator 2015-06-29T10:32:10
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = DensToolKit
TEMPLATE = app


SOURCES += main.cpp\
        dtkmainwindow.cpp \
    dtkglwidget.cpp

HEADERS  += dtkmainwindow.h \
    dtkglwidget.h

FORMS    += dtkmainwindow.ui

unix:!macx {
   LIBS += -lglut
}
macx:  {
   LIBS += -framework GLUT
}
