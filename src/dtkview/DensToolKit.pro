#-------------------------------------------------
#
# Project created by QtCreator 2015-06-29T10:32:10
#
#-------------------------------------------------

QT       += core gui opengl

#greaterThan(QT_MAJOR_VERSION, 5): QT += openglwidgets
greaterThan(QT_MAJOR_VERSION, 5): QT += widgets

TARGET = DensToolKit
TEMPLATE = app


SOURCES += main.cpp \
   dtkmainwindow.cpp \
   dtkglwidget.cpp \
   dtkglbondnetwork.cpp \
   dtkglutils.cpp \
   dtkglcriticalpointnetwork.cpp \
   ../common/bondnetwork.cpp \
   ../common/iofuncts-cpx.cpp \
   ../common/iofuncts-wfn.cpp \
   ../common/iofuncts-wfx.cpp \
   ../common/povraytools.cpp \
   ../common/mymemory.cpp \
   ../common/atomradiicust.cpp \
   ../common/figname.cpp \
   ../common/screenutils.cpp \
   ../common/atomcolschjmol.cpp \
   ../common/stringtools.cpp \
   ../common/critptnetwork.cpp \
   ../common/gausswavefunction.cpp \
   ../common/eigendecompositionjama.cpp \
   ../common/mymath.cpp \
   ../common/fileutils.cpp \
   ../common/matrixvectoroperations3d.cpp \
   ../common/atom.cpp \
   ../common/molecule.cpp \
   ../common/inputmolecule_xyz.cpp

HEADERS  += dtkmainwindow.h \
   dtkglwidget.h \
   localdefs.h \
   dtkglbondnetwork.h \
   dtkglutils.h \
   dtkglcriticalpointnetwork.h \
   ../common/bondnetwork.h \
   ../common/iofuncts-cpx.h \
   ../common/iofuncts-wfn.h \
   ../common/iofuncts-wfx.h \
   ../common/povraytools.h \
   ../common/mymemory.h \
   ../common/atomradiicust.h \
   ../common/screenutils.h \
   ../common/atomcolschjmol.h \
   ../common/stringtools.h \
   ../common/critptnetwork.h \
   ../common/gausswavefunction.h \
   ../common/eigendecompositionjama.h \
   ../common/mymath.h \
   ../common/fileutils.h \
   ../common/figname.h \
   ../common/matrixvectoroperations3d.h \
   ../common/atom.h \
   ../common/molecule.h \
   ../common/inputmolecule_xyz.h

FORMS    += dtkmainwindow.ui

unix:!macx {
   LIBS += -lglut
}
macx:  {
   LIBS += -framework GLUT
}

QMAKE_CXXFLAGS += -include ~/Documents/LongRun/proj/2024dtk/src/dtkview/localdefs.h -include ~/Documents/LongRun/proj/2024dtk/src/common/globaldefs.h -O2

CONFIG += c++11

DISTFILES +=

RESOURCES += \
    denstoolkitviewer.qrc
