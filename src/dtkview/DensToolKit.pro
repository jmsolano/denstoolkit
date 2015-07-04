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
    dtkglwidget.cpp \
    ../common/bondnetwork.cpp \
    ../common/iofuncts-cpx.cpp \
    ../common/iofuncts-wfn.cpp \
    ../common/iofuncts-wfx.cpp \
    ../common/solpovtools.cpp \
    ../common/solmemhand.cpp \
    ../common/atomradiicust.cpp \
    ../common/solscrutils.cpp \
    ../common/atomcolschjmol.cpp \
    ../common/solstringtools.cpp \
    ../common/critptnetwork.cpp \
    ../common/wavefunctionclass.cpp \
    ../common/eig2-4.cpp \
    ../common/solmath.cpp \
    ../common/solfileutils.cpp \
    dtkglbondnetwork.cpp \
    dtkglutils.cpp

HEADERS  += dtkmainwindow.h \
    dtkglwidget.h \
    ../common/bondnetwork.h \
    ../common/iofuncts-cpx.h \
    ../common/iofuncts-wfn.h \
    ../common/iofuncts-wfx.h \
    soldefines.h \
    ../common/solpovtools.h \
    ../common/solmemhand.h \
    ../common/atomradiicust.h \
    ../common/solscrutils.h \
    ../common/atomcolschjmol.h \
    ../common/solstringtools.h \
    ../common/critptnetwork.h \
    ../common/wavefunctionclass.h \
    ../common/eig2-4.h \
    ../common/solmath.h \
    ../common/solfileutils.h \
    dtkglbondnetwork.h \
    dtkglutils.h

FORMS    += dtkmainwindow.ui

unix:!macx {
   LIBS += -lglut
}
macx:  {
   LIBS += -framework GLUT
}

QMAKE_CXXFLAGS += -include soldefines.h
