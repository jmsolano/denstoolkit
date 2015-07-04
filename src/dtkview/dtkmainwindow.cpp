#include "dtkmainwindow.h"
#include "ui_dtkmainwindow.h"
#include <QString>

DTKMainWindow::DTKMainWindow(QWidget *parent) :
   QMainWindow(parent),
   ui(new Ui::DTKMainWindow)
{
   ui->setupUi(this);

   connect(ui->openGLWidget,SIGNAL(rotationChanged()),this,SLOT(updateStatusBar()));
   connect(ui->openGLWidget,SIGNAL(zoomChanged()),this,SLOT(updateStatusBar()));
   updateStatusBar();

   createActions();
   createMenus();
}

void DTKMainWindow::updateStatusBar()
{
   statusBar()->clearMessage();
   statusBar()->showMessage(QString("xAngle: %1   yAngle: %2   zAngle: %3   "
                                    "zoom: %4")
                            .arg(ui->openGLWidget->getXRot())
                            .arg(ui->openGLWidget->getYRot())
                            .arg(ui->openGLWidget->getZRot())
                            .arg(ui->openGLWidget->getCurrentZoom())
                            );
}

DTKMainWindow::~DTKMainWindow()
{
   delete ui;
   delete loadMoleculeAction;
}

void DTKMainWindow::on_resetPushButton_clicked()
{
    ui->openGLWidget->resetView();
}

void DTKMainWindow::createMenus()
{
   fileMenu = menuBar()->addMenu(tr("&File"));
   fileMenu->addAction(loadMoleculeAction);
}

void DTKMainWindow::createActions()
{
    loadMoleculeAction = new QAction(tr("&Open Molecule..."), this);
    loadMoleculeAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_O));
    connect(loadMoleculeAction,SIGNAL(triggered()),this,SLOT(loadMolecule()));
}

void DTKMainWindow::loadMolecule()
{
   QString fname=tr("/Users/jmsolano/Documents/LongRun/proj/readwfn/wavefiles/cubano_sto3g.wfx");
   ui->openGLWidget->addMolecule(fname);
   ui->openGLWidget->update();
}
