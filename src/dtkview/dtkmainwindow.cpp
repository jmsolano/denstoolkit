#include "dtkmainwindow.h"
#include "ui_dtkmainwindow.h"

DTKMainWindow::DTKMainWindow(QWidget *parent) :
   QMainWindow(parent),
   ui(new Ui::DTKMainWindow)
{
   ui->setupUi(this);

   connect(ui->openGLWidget,SIGNAL(rotationChanged()),this,SLOT(updateStatusBar()));
   connect(ui->openGLWidget,SIGNAL(zoomChanged()),this,SLOT(updateStatusBar()));
   updateStatusBar();
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
}

void DTKMainWindow::on_resetPushButton_clicked()
{
    ui->openGLWidget->resetView();
}
