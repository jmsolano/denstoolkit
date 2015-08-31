#include "dtkmainwindow.h"
#include "ui_dtkmainwindow.h"
#include <QString>
#include <QFileDialog>
#include <QImage>
#include <QIcon>
#include <QMessageBox>

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
   setupMainToolbar();
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
   delete loadMoleculeAction;
   delete loadTestMoleculeAction;
   delete clearViewPortAction;
   delete exportViewPortImageAction;
   delete viewAtomLabelsAction;
   delete showAboutDTKAction;

   delete fileMenu;
   delete viewMenu;
   delete helpMenu;

   delete fileToolBar;

   delete ui;
}

void DTKMainWindow::on_resetPushButton_clicked()
{
    ui->openGLWidget->resetView();
}

void DTKMainWindow::createMenus()
{
#ifdef __linux__
   menuBar()->setNativeMenuBar(false);
#endif
   fileMenu = menuBar()->addMenu(tr("&File"));
   fileMenu->addAction(loadMoleculeAction);
   fileMenu->addAction(loadTestMoleculeAction);
   fileMenu->addAction(clearViewPortAction);
   fileMenu->addAction(exportViewPortImageAction);

   viewMenu = menuBar()->addMenu(tr("&View"));
   viewMenu->addAction(viewAtomLabelsAction);

   helpMenu = menuBar()->addMenu(tr("&Help"));
   helpMenu->addAction(showAboutDTKAction);
}

void DTKMainWindow::createActions()
{
    loadMoleculeAction = new QAction(QIcon(":/images/open.png"),tr("&Open Molecule..."), this);
    loadMoleculeAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_O));
    connect(loadMoleculeAction,SIGNAL(triggered()),this,SLOT(loadMolecule()));

    loadTestMoleculeAction = new QAction(tr("Open &Test Molecule"), this);
    loadTestMoleculeAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_T));
    connect(loadTestMoleculeAction,SIGNAL(triggered()),this,SLOT(loadTestMolecule()));

    clearViewPortAction = new QAction(QIcon(":/images/square.png"),tr("&Clear Molecules"), this);
    clearViewPortAction->setShortcut(QKeySequence(Qt::CTRL+Qt::ShiftModifier+Qt::Key_C));
    connect(clearViewPortAction,SIGNAL(triggered()),this,SLOT(clearViewPort()));

    exportViewPortImageAction = new QAction(QIcon(":/images/save.png"),tr("&Export Image..."), this);
    exportViewPortImageAction->setShortcut(QKeySequence(Qt::CTRL+Qt::ShiftModifier+Qt::Key_E));
    connect(exportViewPortImageAction,SIGNAL(triggered()),this,SLOT(exportViewPortImage()));

    viewAtomLabelsAction = new QAction(QIcon(":/images/showatlbls.png"),tr("View Atom &Labels"), this);
    viewAtomLabelsAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_L));
    viewAtomLabelsAction->setCheckable(true);
    viewAtomLabelsAction->setChecked(ui->viewAtLblsCheckBox->isChecked());
    connect(viewAtomLabelsAction,SIGNAL(triggered()),this,SLOT(setViewAtomLabels()));

    showAboutDTKAction = new QAction(tr("&About DTK"), this);
    connect(showAboutDTKAction,SIGNAL(triggered()),this,SLOT(showAboutDTK()));
}

void DTKMainWindow::setupMainToolbar()
{
   ui->mainToolBar->addAction(loadMoleculeAction);
   ui->mainToolBar->addAction(exportViewPortImageAction);
   ui->mainToolBar->addAction(clearViewPortAction);
   ui->mainToolBar->addSeparator();
   ui->mainToolBar->addAction(viewAtomLabelsAction);
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

void DTKMainWindow::setViewAtomLabels()
{
   bool val=ui->viewAtLblsCheckBox->isChecked();
   val=(!val);
   ui->viewAtLblsCheckBox->setChecked(val);
   viewAtomLabelsAction->setChecked(val);
   ui->openGLWidget->setDrawAtomLabels(val);
}

void DTKMainWindow::showAboutDTK()
{
   QMessageBox about(this);
   about.setText(tr("DensToolKitViewer"));
   about.setInformativeText("(C) J. M. Solano-Altamirano, 2015."
                            "\nThis program is licensed with a GPLv3"
                            "\nlicense. See\n\nhttp://www.gnu.org/licenses/gpl-3.0.en.html\n"
                            "\nfor more information.");
   QPixmap myicon(QString(":/images/dtk256x256.png"));
   about.setIconPixmap(myicon);
   int ret=about.exec();
}

void DTKMainWindow::on_viewAtLblsCheckBox_clicked()
{
   bool val=ui->viewAtLblsCheckBox->isChecked();
    ui->openGLWidget->setDrawAtomLabels(val);
    viewAtomLabelsAction->setChecked(val);
}
