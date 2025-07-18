/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 2.1.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2025, Juan Manuel Solano-Altamirano
                                   <jmsolanoalt@gmail.com>
  
   -------------------------------------------------------------------
  
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
   ---------------------------------------------------------------------
  
   If you want to redistribute modifications of the suite, please
   consider to include your modifications in our official release.
   We will be pleased to consider the inclusion of your code
   within the official distribution. Please keep in mind that
   scientific software is very special, and version control is 
   crucial for tracing bugs. If in despite of this you distribute
   your modified version, please do not call it DensToolKit.
  
   If you find DensToolKit useful, we humbly ask that you cite
   the paper(s) on the package --- you can find them on the top
   README file.
*/

#include "dtkmainwindow.h"
#include "ui_dtkmainwindow.h"
#include <QString>
#include <QFileDialog>
#include <QImage>
#include <QIcon>
#include <QMessageBox>

DTKMainWindow::DTKMainWindow(QWidget *parent) :
   QMainWindow(parent),
   ui(new Ui::DTKMainWindow) {
   ui->setupUi(this);

   connect(ui->openGLWidget,SIGNAL(rotationChanged()),this,SLOT(updateStatusBar()));
   connect(ui->openGLWidget,SIGNAL(rotationChanged()),this,SLOT(updateAngleScrollBarsAndSpinBoxes()));
   connect(ui->openGLWidget,SIGNAL(zoomChanged()),this,SLOT(updateStatusBar()));
   updateStatusBar();

   createActions();
   createMenus();
   setupMainToolbar();
   ui->menuBar->setNativeMenuBar(true);
}
void DTKMainWindow::updateStatusBar() {
   statusBar()->clearMessage();
   statusBar()->showMessage(QString("xAngle: %1   yAngle: %2   zAngle: %3   "
                                    "zoom: %4")
                            .arg(ui->openGLWidget->getXRot())
                            .arg(ui->openGLWidget->getYRot())
                            .arg(ui->openGLWidget->getZRot())
                            .arg(ui->openGLWidget->getCurrentZoom())
                            );
}
void DTKMainWindow::updateAngleScrollBarsAndSpinBoxes() {
   ui->xAngleSpinBox->setValue(ui->openGLWidget->getXRot());
   ui->yAngleSpinBox->setValue(ui->openGLWidget->getYRot());
   ui->zAngleSpinBox->setValue(ui->openGLWidget->getZRot());

   ui->xAngleScrollBar->setValue(ui->openGLWidget->getXRot());
   ui->yAngleScrollBar->setValue(ui->openGLWidget->getYRot());
   ui->zAngleScrollBar->setValue(ui->openGLWidget->getZRot());
}
DTKMainWindow::~DTKMainWindow() {
   delete loadMoleculeAction;
   delete loadTestMoleculeAction;
   delete clearViewPortAction;
   delete exportViewPortImageAction;

   delete viewAtomsAction;
   delete viewAtomLabelsAction;
   delete viewRegularBondsAction;
   delete setTransparentAtomsAndLinksAction;

   delete viewBondGradientPathsAction;
   delete viewRingGradientPathsAction;
   delete viewCageGradientPathsAction;
   delete viewACPLabelsAction;
   delete viewBCPLabelsAction;
   delete viewRCPLabelsAction;
   delete viewCCPLabelsAction;

   delete showAboutDTKAction;

   delete fileMenu;
   delete viewMenu;
   delete helpMenu;

   //delete fileToolBar;

   delete ui;
}
void DTKMainWindow::on_resetPushButton_clicked() {
    ui->openGLWidget->resetView();
}
void DTKMainWindow::createMenus() {
#ifdef __linux__
   menuBar()->setNativeMenuBar(false);
#endif
   fileMenu = menuBar()->addMenu(tr("&File"));
   fileMenu->addAction(loadMoleculeAction);
   fileMenu->addAction(loadTestMoleculeAction);
   fileMenu->addAction(clearViewPortAction);
   fileMenu->addAction(exportViewPortImageAction);

   viewMenu = menuBar()->addMenu(tr("&View options"));
   viewMenu->addAction(viewAtomsAction);
   viewMenu->addAction(viewRegularBondsAction);
   viewMenu->addAction(setTransparentAtomsAndLinksAction);
   viewMenu->addAction(viewAtomLabelsAction);
   viewMenu->addSeparator();
   viewMenu->addAction(viewBondGradientPathsAction);
   viewMenu->addAction(viewRingGradientPathsAction);
   viewMenu->addAction(viewCageGradientPathsAction);
   viewMenu->addAction(viewACPLabelsAction);
   viewMenu->addAction(viewBCPLabelsAction);
   viewMenu->addAction(viewRCPLabelsAction);
   viewMenu->addAction(viewCCPLabelsAction);

   helpMenu = menuBar()->addMenu(tr("&Help"));
   helpMenu->addAction(showAboutDTKAction);
}
void DTKMainWindow::createActions() {
    loadMoleculeAction = new QAction(QIcon(":/images/open.png"),tr("&Open Molecule..."), this);
    loadMoleculeAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_O));
    connect(loadMoleculeAction,SIGNAL(triggered()),this,SLOT(loadMolecule()));

    loadTestMoleculeAction = new QAction(tr("Open &Test Molecule"), this);
    loadTestMoleculeAction->setShortcut(QKeySequence(Qt::CTRL+Qt::ShiftModifier+Qt::Key_T));
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

    viewBondGradientPathsAction = new QAction(QIcon(":/images/viewbgps.png"),tr("Draw &Bond Gradient Paths"), this);
    viewBondGradientPathsAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_1));
    viewBondGradientPathsAction->setCheckable(true);
    viewBondGradientPathsAction->setChecked(ui->viewBGPsCheckBox->isChecked());
    connect(viewBondGradientPathsAction,SIGNAL(triggered()),this,SLOT(setViewBondGradientPaths()));

    viewRingGradientPathsAction = new QAction(QIcon(":/images/viewrgps.png"),tr("Draw &Ring Gradient Paths"), this);
    viewRingGradientPathsAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_2));
    viewRingGradientPathsAction->setCheckable(true);
    viewRingGradientPathsAction->setChecked(ui->viewRGPsCheckBox->isChecked());
    connect(viewRingGradientPathsAction,SIGNAL(triggered()),this,SLOT(setViewRingGradientPaths()));

    viewCageGradientPathsAction = new QAction(QIcon(":/images/viewcgps.png"),tr("Draw &Cage Gradient Paths"), this);
    viewCageGradientPathsAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_3));
    viewCageGradientPathsAction->setCheckable(true);
    viewCageGradientPathsAction->setChecked(ui->viewCGPsCheckBox->isChecked());
    connect(viewCageGradientPathsAction,SIGNAL(triggered()),this,SLOT(setViewCageGradientPaths()));

    //viewACPLabelsAction = new QAction(QIcon(":/images/drawcplbls.png"),tr("View C&P Labels"), this);
    viewACPLabelsAction = new QAction(QIcon(":/images/drawacplbls.png"),tr("View ACP Labels"), this);
    //viewACPLabelsAction->setShortcut(QKeySequence(Qt::CTRL+Qt::ShiftModifier+Qt::Key_L));
    viewACPLabelsAction->setCheckable(true);
    viewACPLabelsAction->setChecked(ui->viewACPLabelsCheckBox->isChecked());
    connect(viewACPLabelsAction,SIGNAL(triggered()),this,SLOT(setViewACPLabels()));

    viewBCPLabelsAction = new QAction(QIcon(":/images/drawbcplbls.png"),tr("View BCP Labels"), this);
    //viewBCPLabelsAction->setShortcut(QKeySequence(Qt::CTRL+Qt::ShiftModifier+Qt::Key_L));
    viewBCPLabelsAction->setCheckable(true);
    viewBCPLabelsAction->setChecked(ui->viewBCPLabelsCheckBox->isChecked());
    connect(viewBCPLabelsAction,SIGNAL(triggered()),this,SLOT(setViewBCPLabels()));

    viewRCPLabelsAction = new QAction(QIcon(":/images/drawrcplbls.png"),tr("View RCP Labels"), this);
    //viewRCPLabelsAction->setShortcut(QKeySequence(Qt::CTRL+Qt::ShiftModifier+Qt::Key_L));
    viewRCPLabelsAction->setCheckable(true);
    viewRCPLabelsAction->setChecked(ui->viewRCPLabelsCheckBox->isChecked());
    connect(viewRCPLabelsAction,SIGNAL(triggered()),this,SLOT(setViewRCPLabels()));

    viewCCPLabelsAction = new QAction(QIcon(":/images/drawccplbls.png"),tr("View CCP Labels"), this);
    //viewCCPLabelsAction->setShortcut(QKeySequence(Qt::CTRL+Qt::ShiftModifier+Qt::Key_L));
    viewCCPLabelsAction->setCheckable(true);
    viewCCPLabelsAction->setChecked(ui->viewCCPLabelsCheckBox->isChecked());
    connect(viewCCPLabelsAction,SIGNAL(triggered()),this,SLOT(setViewCCPLabels()));

    setTransparentAtomsAndLinksAction = new QAction(QIcon(":/images/toggletransp.png"),tr("Set &Transparent"), this);
    setTransparentAtomsAndLinksAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_T));
    setTransparentAtomsAndLinksAction->setCheckable(true);
    setTransparentAtomsAndLinksAction->setChecked(ui->setTransparentCheckBox->isChecked());
    connect(setTransparentAtomsAndLinksAction,SIGNAL(triggered()),this,SLOT(setTransparentAtomsAndLinks()));

    viewRegularBondsAction = new QAction(QIcon(":/images/drawlinks.png"),tr("View R&egular Bonds"), this);
    viewRegularBondsAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_B));
    viewRegularBondsAction->setCheckable(true);
    viewRegularBondsAction->setChecked(ui->viewAtLblsCheckBox->isChecked());
    connect(viewRegularBondsAction,SIGNAL(triggered()),this,SLOT(setViewRegularBonds()));

    viewAtomsAction = new QAction(QIcon(":/images/drawatoms.png"),tr("View &Atoms"), this);
    viewAtomsAction->setShortcut(QKeySequence(Qt::CTRL+Qt::Key_A));
    viewAtomsAction->setCheckable(true);
    viewAtomsAction->setChecked(ui->viewAtomsCheckBox->isChecked());
    connect(viewAtomsAction,SIGNAL(triggered()),this,SLOT(setViewAtoms()));

    showAboutDTKAction = new QAction(tr("&About DTK"), this);
    connect(showAboutDTKAction,SIGNAL(triggered()),this,SLOT(showAboutDTK()));
}
void DTKMainWindow::setupMainToolbar() {
   ui->mainToolBar->addAction(loadMoleculeAction);
   ui->mainToolBar->addAction(exportViewPortImageAction);
   ui->mainToolBar->addAction(clearViewPortAction);
   ui->mainToolBar->addSeparator();
   ui->mainToolBar->addAction(viewAtomsAction);
   ui->mainToolBar->addAction(viewRegularBondsAction);
   ui->mainToolBar->addAction(viewAtomLabelsAction);
   ui->mainToolBar->addAction(setTransparentAtomsAndLinksAction);
   ui->mainToolBar->addSeparator();
   ui->mainToolBar->addAction(viewBondGradientPathsAction);
   ui->mainToolBar->addAction(viewRingGradientPathsAction);
   ui->mainToolBar->addAction(viewCageGradientPathsAction);
   ui->mainToolBar->addAction(viewACPLabelsAction);
   ui->mainToolBar->addAction(viewBCPLabelsAction);
   ui->mainToolBar->addAction(viewRCPLabelsAction);
   ui->mainToolBar->addAction(viewCCPLabelsAction);
}
void DTKMainWindow::loadMolecule() {
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
void DTKMainWindow::loadTestMolecule() {
#ifdef __APPLE__
   QString fname=tr("/Users/jmsolano/Documents/LongRun/proj/readwfn/src/dtkview/cubano_sto3gRhoCP.cpx");
   //QString fname=tr("/Users/jmsolano/Documents/LongRun/proj/readwfn/src/dtkview/phenantreneRhoCP.cpx");
#else
   QString fname=tr("/home/jmsolano/Documents/prog/dtk/src/dtkview/phenantreneRhoCP.cpx");
   //QString fname=tr("/home/jmsolano/Documents/prog/dtk/src/dtkview/cubano_sto3gRhoCP.cpx");
#endif
   ui->openGLWidget->addMolecule(fname);
   ui->openGLWidget->update();
}
void DTKMainWindow::clearViewPort() {
   ui->openGLWidget->clearWFsBNsCPXs();
   ui->openGLWidget->update();
}
void DTKMainWindow::exportViewPortImage() {
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
void DTKMainWindow::setTransparentAtomsAndLinks() {
   bool val=ui->setTransparentCheckBox->isChecked();
   val=(!val);
   ui->setTransparentCheckBox->setChecked(val);
   setTransparentAtomsAndLinksAction->setChecked(val);
   ui->openGLWidget->setTransparentAtomsAndLinks(val);
}
void DTKMainWindow::setViewAtoms() {
   bool val=ui->viewAtomsCheckBox->isChecked();
   val=(!val);
   ui->viewAtomsCheckBox->setChecked(val);
   viewAtomsAction->setChecked(val);
   ui->openGLWidget->setViewAtoms(val);
}
void DTKMainWindow::setViewAtomLabels() {
   bool val=ui->viewAtLblsCheckBox->isChecked();
   val=(!val);
   ui->viewAtLblsCheckBox->setChecked(val);
   viewAtomLabelsAction->setChecked(val);
   ui->openGLWidget->setDrawAtomLabels(val);
}
void DTKMainWindow::setViewRegularBonds() {
   bool val=ui->viewRegularBondsCheckBox->isChecked();
   val=(!val);
   ui->viewRegularBondsCheckBox->setChecked(val);
   viewRegularBondsAction->setChecked(val);
   ui->openGLWidget->setViewRegularBonds(val);
}
void DTKMainWindow::setViewBondGradientPaths() {
   bool val=ui->viewBGPsCheckBox->isChecked();
   val=(!val);
   ui->viewBGPsCheckBox->setChecked(val);
   viewBondGradientPathsAction->setChecked(val);
   ui->openGLWidget->setViewBondGradientPaths(val);
}
void DTKMainWindow::setViewRingGradientPaths() {
   bool val=ui->viewRGPsCheckBox->isChecked();
   val=(!val);
   ui->viewRGPsCheckBox->setChecked(val);
   viewRingGradientPathsAction->setChecked(val);
   ui->openGLWidget->setViewRingGradientPaths(val);
}
void DTKMainWindow::setViewCageGradientPaths() {
   bool val=ui->viewCGPsCheckBox->isChecked();
   val=(!val);
   ui->viewCGPsCheckBox->setChecked(val);
   viewCageGradientPathsAction->setChecked(val);
   ui->openGLWidget->setViewCageGradientPaths(val);
}
void DTKMainWindow::setViewACPLabels() {
   bool val=ui->viewACPLabelsCheckBox->isChecked();
   val=(!val);
   ui->viewACPLabelsCheckBox->setChecked(val);
   viewACPLabelsAction->setChecked(val);
   ui->openGLWidget->setDrawACPLabels(val);
}
void DTKMainWindow::setViewBCPLabels() {
   bool val=ui->viewBCPLabelsCheckBox->isChecked();
   val=(!val);
   ui->viewBCPLabelsCheckBox->setChecked(val);
   viewBCPLabelsAction->setChecked(val);
   ui->openGLWidget->setDrawBCPLabels(val);
}
void DTKMainWindow::setViewRCPLabels() {
   bool val=ui->viewRCPLabelsCheckBox->isChecked();
   val=(!val);
   ui->viewRCPLabelsCheckBox->setChecked(val);
   viewRCPLabelsAction->setChecked(val);
   ui->openGLWidget->setDrawRCPLabels(val);
}
void DTKMainWindow::setViewCCPLabels() {
   bool val=ui->viewCCPLabelsCheckBox->isChecked();
   val=(!val);
   ui->viewCCPLabelsCheckBox->setChecked(val);
   viewCCPLabelsAction->setChecked(val);
   ui->openGLWidget->setDrawCCPLabels(val);
}
void DTKMainWindow::showAboutDTK() {
   QMessageBox about(NULL);
   about.setText(tr("About DensToolKitViewer"));
   about.setInformativeText(tr("Version: "
                            CURRENTVERSION 
                            "\n\n(C) J. M. Solano-Altamirano, 2015-2020."
                            "\nThis program is licensed with a GPLv3"
                            "\nlicense. See\n\nhttp://www.gnu.org/licenses/gpl-3.0.en.html\n"
                            "\nfor more information."
                            "\n\nIf you use this program, please consider citing our paper,\n"
                            "we will appreciate it very much.\n\n"
                            "Comput. Phys. Commun., doi: 10.1016/j.cpc.2015.07.005"));
   QPixmap myicon(QString(":/images/dtk256x256.png"));
   about.setIconPixmap(myicon);
   int ret=about.exec();
}
void DTKMainWindow::on_viewAtLblsCheckBox_clicked() {
   bool val=ui->viewAtLblsCheckBox->isChecked();
    ui->openGLWidget->setDrawAtomLabels(val);
    viewAtomLabelsAction->setChecked(val);
}
void DTKMainWindow::on_viewCGPsCheckBox_clicked() {
   bool val=ui->viewCGPsCheckBox->isChecked();
   ui->openGLWidget->setViewCageGradientPaths(val);
   viewBondGradientPathsAction->setChecked(val);
}
void DTKMainWindow::on_viewRGPsCheckBox_clicked() {
   bool val=ui->viewRGPsCheckBox->isChecked();
   ui->openGLWidget->setViewRingGradientPaths(val);
   viewRingGradientPathsAction->setChecked(val);
}
void DTKMainWindow::on_viewBGPsCheckBox_clicked() {
   bool val=ui->viewBGPsCheckBox->isChecked();
   ui->openGLWidget->setViewBondGradientPaths(val);
   viewCageGradientPathsAction->setChecked(val);
}
void DTKMainWindow::on_viewRegularBondsCheckBox_clicked() {
   bool val=ui->viewRegularBondsCheckBox->isChecked();
   ui->openGLWidget->setViewRegularBonds(val);
   viewRegularBondsAction->setChecked(val);
}
void DTKMainWindow::on_viewAtomsCheckBox_clicked() {
   bool val=ui->viewAtomsCheckBox->isChecked();
   ui->openGLWidget->setViewAtoms(val);
   viewAtomsAction->setChecked(val);
}
void DTKMainWindow::on_setTransparentCheckBox_clicked() {
   bool val=ui->setTransparentCheckBox->isChecked();
   ui->openGLWidget->setTransparentAtomsAndLinks(val);
   setTransparentAtomsAndLinksAction->setChecked(val);
}
void DTKMainWindow::on_viewACPLabelsCheckBox_clicked() {
   bool val=ui->viewACPLabelsCheckBox->isChecked();
   ui->openGLWidget->setDrawACPLabels(val);
   viewACPLabelsAction->setChecked(val);
}
void DTKMainWindow::on_viewBCPLabelsCheckBox_clicked() {
   bool val=ui->viewBCPLabelsCheckBox->isChecked();
   ui->openGLWidget->setDrawBCPLabels(val);
   viewBCPLabelsAction->setChecked(val);
}
void DTKMainWindow::on_viewRCPLabelsCheckBox_clicked() {
   bool val=ui->viewRCPLabelsCheckBox->isChecked();
   ui->openGLWidget->setDrawRCPLabels(val);
   viewRCPLabelsAction->setChecked(val);
}
void DTKMainWindow::on_viewCCPLabelsCheckBox_clicked() {
   bool val=ui->viewCCPLabelsCheckBox->isChecked();
   ui->openGLWidget->setDrawCCPLabels(val);
   viewCCPLabelsAction->setChecked(val);
}
void DTKMainWindow::on_xAngleScrollBar_valueChanged(int value) {
   ui->xAngleSpinBox->setValue(value);
   ui->openGLWidget->setXRotation(value);
}
void DTKMainWindow::on_yAngleScrollBar_valueChanged(int value) {
   ui->yAngleSpinBox->setValue(value);
   ui->openGLWidget->setYRotation(value);
}
void DTKMainWindow::on_zAngleScrollBar_valueChanged(int value) {
   ui->zAngleSpinBox->setValue(value);
   ui->openGLWidget->setZRotation(value);
}
void DTKMainWindow::on_xAngleSpinBox_valueChanged(int value) {
   ui->xAngleScrollBar->setValue(value);
   ui->openGLWidget->setXRotation(value);
}
void DTKMainWindow::on_yAngleSpinBox_valueChanged(int value) {
   ui->yAngleScrollBar->setValue(value);
   ui->openGLWidget->setYRotation(value);
}
void DTKMainWindow::on_zAngleSpinBox_valueChanged(int value) {
   ui->zAngleScrollBar->setValue(value);
   ui->openGLWidget->setZRotation(value);
}

