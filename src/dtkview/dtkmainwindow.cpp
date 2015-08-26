#include "dtkmainwindow.h"
#include "ui_dtkmainwindow.h"
#include <QString>
#include <QFileDialog>
#include <QImage>

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
   fileMenu->addAction(loadTestMoleculeAction);
   fileMenu->addAction(clearViewPortAction);
   fileMenu->addAction(exportViewPortImageAction);
}

void DTKMainWindow::createActions()
{
    loadMoleculeAction = new QAction(tr("&Open Molecule..."), this);
    loadMoleculeAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_O));
    connect(loadMoleculeAction,SIGNAL(triggered()),this,SLOT(loadMolecule()));

    loadTestMoleculeAction = new QAction(tr("Open &Test Molecule"), this);
    loadTestMoleculeAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_T));
    connect(loadTestMoleculeAction,SIGNAL(triggered()),this,SLOT(loadTestMolecule()));

    clearViewPortAction = new QAction(tr("&Clear Molecules"), this);
    clearViewPortAction->setShortcut(QKeySequence(Qt::CTRL+Qt::ShiftModifier+Qt::Key_C));
    connect(clearViewPortAction,SIGNAL(triggered()),this,SLOT(clearViewPort()));

    exportViewPortImageAction = new QAction(tr("&Export Image..."), this);
    exportViewPortImageAction->setShortcut(QKeySequence(Qt::CTRL+Qt::ShiftModifier+Qt::Key_E));
    connect(exportViewPortImageAction,SIGNAL(triggered()),this,SLOT(exportViewPortImage()));
}

void DTKMainWindow::loadMolecule()
{
#ifdef __APPLE__
   QString fname=tr("/Users/jmsolano/Documents/LongRun/proj/readwfn/src/dtkview/cubano_sto3gRhoCP.cpx");
   //QString fname=tr("/Users/jmsolano/Documents/LongRun/proj/readwfn/src/dtkview/phenantreneRhoCP.cpx");
#else
   //QString fname=tr("/home/jmsolano/Documents/prog/dtk/src/dtkview/phenantreneRhoCP.cpx");
   QString fname=tr("/home/jmsolano/Documents/prog/dtk/src/dtkview/cubano_sto3gRhoCP.cpx");
#endif
   QString selfilter = tr("DTK CPX files (*.cpx)");
   QFileDialog *myfdiag=new QFileDialog;
   fname=myfdiag->getOpenFileName(this,
                                      tr("Load cpx file"),
                                      QDir::homePath(),
                                      tr("DTK CPX files (*.cpx)"\
                                      /*";;WFN files (*.wfn)"\
                                      ";;WFX files (*.wfx)"\
                                         */),
                                      &selfilter
                                      );
   //qDebug() << "fname: " << fname;
   delete myfdiag;
   if (fname.size()==0) { return; }
   ui->openGLWidget->addMolecule(fname);
   ui->openGLWidget->update();
}

void DTKMainWindow::loadTestMolecule()
{
#ifdef __APPLE__
   //QString fname=tr("/Users/jmsolano/Documents/LongRun/proj/readwfn/src/dtkview/cubano_sto3gRhoCP.cpx");
   QString fname=tr("/Users/jmsolano/Documents/LongRun/proj/readwfn/src/dtkview/phenantreneRhoCP.cpx");
#else
   QString fname=tr("/home/jmsolano/Documents/prog/dtk/src/dtkview/phenantreneRhoCP.cpx");
   //QString fname=tr("/home/jmsolano/Documents/prog/dtk/src/dtkview/cubano_sto3gRhoCP.cpx");
#endif
   ui->openGLWidget->addMolecule(fname);
   ui->openGLWidget->update();
}

void DTKMainWindow::clearViewPort()
{
   ui->openGLWidget->clearWFsBNsCPXs();
   ui->openGLWidget->update();
}

void DTKMainWindow::exportViewPortImage()
{
    QString selfilter = tr("Portable Network Graphics (*.png)");
    QFileDialog *myfdiag=new QFileDialog;
    QString fname=myfdiag->getSaveFileName(this,
                                       tr("Save Image"),
                                       QDir::homePath(),
                                       tr("Portable Network Graphcis (*.png)"\
                                       ";;JPEG (*.jpg)"\
                                       ";;Windows BitMap (*.bmp)"\
                                       ";;All suported images (*.bmp *.png *.jpg)"
                                          ),
                                       &selfilter
                                       );
    delete myfdiag;
    if (fname.size()==0) { return; }
#ifdef __APPLE__
    QImage imgbuff=ui->openGLWidget->grabFramebuffer();
#else
    QImage imgbuff=ui->openGLWidget->grabFrameBuffer(true);
#endif
    imgbuff.save(fname,0,100);
}

void DTKMainWindow::on_viewAtLblsCheckBox_clicked()
{
    ui->openGLWidget->setDrawAtomLabels(ui->viewAtLblsCheckBox->isChecked());
}
