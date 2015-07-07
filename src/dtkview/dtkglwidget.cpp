#include "dtkglwidget.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <QOpenGLWidget>
#else
#include <GL/glut.h>
#include <QGLWidget>
#endif
#include <QtGui>
#include <QMessageBox>
#include <QVector3D>
#include <QMatrix4x4>
#include <QDebug>
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
#include "dtkglutils.h"
#include "../common/iofuncts-cpx.h"

#define DTKGL_DEFAULT_SPHERE_RADIUS 0.5f
#define DTKGL_DEFAULT_LINK_RADIUS 0.2f
#define DTKGL_DEFAULT_CRIT_PT_RADIUS 0.125f
#define DTKGL_DEFAULT_GRAD_PATH_RADIUS 0.045f

DTKGLWidget::DTKGLWidget(QWidget *parent)
#ifdef __APPLE__
   : QOpenGLWidget(parent)
#else
   : QGLWidget(parent)
   #endif
{
   xRot=0;
   yRot=0;
   zRot=0;
   cameraDistance=INITIAL_CAMERA_DISTANCE;
   waveFunction.clear();
   bondNW.clear();
   critPtNW.clear();
   //zNear=0.01f;
   //zFar=100.0f;
}

DTKGLWidget::~DTKGLWidget()
{
   for (int i=0; i<waveFunction.size(); ++i) {
      delete waveFunction[i];
      waveFunction[i]=NULL;
   }
   for (int i=0; i<bondNW.size(); ++i){
      delete bondNW[i];
      bondNW[i]=NULL;
   }
   for (int i=0; i<critPtNW.size(); ++i) {
      delete critPtNW[i];
      critPtNW[i]=NULL;
   }
}

void DTKGLWidget::initializeGL()
{
   /*
   glClearColor(0.2,0.2,0.2,1.0);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_COLOR_MATERIAL);
   // */

   GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
   GLfloat mat_shininess[] = { 80.0 };
   GLfloat light_position[] = { 2.0f*INITIAL_CAMERA_DISTANCE,\
                                2.0f*INITIAL_CAMERA_DISTANCE,\
                                2.0f*INITIAL_CAMERA_DISTANCE, 0.0 };
   glClearColor((48.f/255.f),(62.f/255.f),(115.f/255.f),1.0f);
   glShadeModel (GL_SMOOTH);

   glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);

   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_COLOR_MATERIAL);
   glEnable(GL_NORMALIZE); //Needed for keeping light uniform
}

void DTKGLWidget::paintGL()
{

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glLoadIdentity();
   glTranslatef(0.0, 0.0, -10.0);
   glRotatef(float(xRot), 1.0, 0.0, 0.0);
   glRotatef(float(yRot), 0.0, 1.0, 0.0);
   glRotatef(float(zRot), 0.0, 0.0, 1.0);
   glScalef(cameraDistance,cameraDistance,cameraDistance);

   drawEverything();
}

void DTKGLWidget::resizeGL(int w, int h)
{
   /*
   glViewport(0,0,w,h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective(40,(float(w)/float(h)),0.001,1000);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   gluLookAt(0,0,cameraDistance, 0,0,0, 0,1,0);
   // */
   double side = qMin(w, h);

   glViewport(0, 0, w, h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glFrustum(-w/side, w/side, -h/side, h/side, 2.5, 200.0);
   glMatrixMode(GL_MODELVIEW);
}

void DTKGLWidget::drawAtoms()
{
   DTKGLBondNetWork *bn;
   QVector3D pos,col;
   float rad = DTKGL_DEFAULT_SPHERE_RADIUS;
   for (int i=0; i<bondNW.size(); ++i) {
      bn=bondNW[i];
      for (int j=0; j<(bn->numAtoms()); ++j) {
         pos=bn->getAtomCoordinates(j);
         //col[0]=0.8f; col[1]=0.0f; col[2]=0.0f;
         col=bn->getAtomColor(j);
         drawSingleSphere(
                     pos[0],pos[1],pos[2],
                     rad, col[0], col[1], col[2]
                 );
      }
   }
}

void DTKGLWidget::drawLinks()
{
    DTKGLBondNetWork *bn;
    float rad = DTKGL_DEFAULT_LINK_RADIUS;
    for (int bIdx=0; bIdx<bondNW.size(); ++bIdx) {
       bn=bondNW[bIdx];
       for (int i=0; i<bn->numLinks(); ++i) {
          drawSingleCylinder(bn->getLinkStart(i),bn->getLinkHeight(i),rad,\
                             bn->getLinkAngle(i),bn->getLinkRotationVector(i),\
                             bn->getLinkColor(i));
       }

    }
}

void DTKGLWidget::drawACPs()
{
   DTKGLCriticalPointNetWork *cp;
   const static QVector3D colACP(0.0f,0.0f,0.0f);
   QVector3D pos;
   float rad=DTKGL_DEFAULT_CRIT_PT_RADIUS;
   for (int cpIdx=0; cpIdx<critPtNW.size(); ++cpIdx) {
      cp=critPtNW[cpIdx];
      if (!(cp->iKnowACPs())) {return;}
      for (int i=0; i<cp->getNumACPs(); ++i) {
         pos=cp->getACPCoordinates(i);
         drawSingleSphere(pos,rad,colACP);
      }
   }
}

void DTKGLWidget::drawBCPs()
{
   DTKGLCriticalPointNetWork *cp;
   const static QVector3D colBCP(0.0f,0.6f,1.0f);
   QVector3D pos;
   float rad=DTKGL_DEFAULT_CRIT_PT_RADIUS;
   for (int cpIdx=0; cpIdx<critPtNW.size(); ++cpIdx) {
      cp=critPtNW[cpIdx];
      if (!(cp->iKnowBCPs())) {return;}
      for (int i=0; i<cp->getNumBCPs(); ++i) {
         pos=cp->getBCPCoordinates(i);
         drawSingleSphere(pos,rad,colBCP);
      }
   }
}

void DTKGLWidget::drawRCPs()
{
   DTKGLCriticalPointNetWork *cp;
   const static QVector3D colRCP(1.0f,1.0f,0.0f);
   QVector3D pos;
   float rad=DTKGL_DEFAULT_CRIT_PT_RADIUS;
   for (int cpIdx=0; cpIdx<critPtNW.size(); ++cpIdx) {
      cp=critPtNW[cpIdx];
      if (!(cp->iKnowRCPs())) {return;}
      for (int i=0; i<cp->getNumRCPs(); ++i) {
         pos=cp->getRCPCoordinates(i);
         drawSingleSphere(pos,rad,colRCP);
      }
   }
}

void DTKGLWidget::drawCCPs()
{
   DTKGLCriticalPointNetWork *cp;
   const static QVector3D colCCP(1.0f,0.0f,0.0f);
   QVector3D pos;
   float rad=DTKGL_DEFAULT_CRIT_PT_RADIUS;
   for (int cpIdx=0; cpIdx<critPtNW.size(); ++cpIdx) {
      cp=critPtNW[cpIdx];
      if (!(cp->iKnowCCPs())) {return;}
      for (int i=0; i<cp->getNumCCPs(); ++i) {
         pos=cp->getCCPCoordinates(i);
         drawSingleSphere(pos,rad,colCCP);
      }
   }
}

void DTKGLWidget::drawCriticalPoints()
{
   drawACPs();
   drawBCPs();
   drawRCPs();
   drawCCPs();
}

void DTKGLWidget::drawBGPs()
{
   DTKGLCriticalPointNetWork *cp;
   const static QVector3D colBGP(0.0f,0.2f,1.0f);
   QVector3D pos;
   int npts;
   float rad=DTKGL_DEFAULT_GRAD_PATH_RADIUS;
   for (int cpIdx=0; cpIdx<critPtNW.size(); ++cpIdx) {
      cp=critPtNW[cpIdx];
      for (int i=0; i<cp->getNumBCPs(); ++i) {
         npts=cp->getNumPtsOfBGP(i);
         for (int j=0; j<npts; ++j) {
            pos=cp->getBGPPointCoordinates(i,j);
            drawSingleSphere(pos,rad,colBGP);
         }
      }
   }
}

void DTKGLWidget::drawRGPs()
{
   DTKGLCriticalPointNetWork *cp;
   const static QVector3D colRGP(0.0f,0.2f,1.0f);
   QVector3D pos;
   int npts,k;
   float rad=DTKGL_DEFAULT_GRAD_PATH_RADIUS;
   for (int cpIdx=0; cpIdx<critPtNW.size(); ++cpIdx) {
      cp=critPtNW[cpIdx];
      for (int i=0; i<cp->getNumRCPs(); ++i) {
         k=0;
      }
   }
}

void DTKGLWidget::drawGradientPaths()
{
   drawBGPs();
}

void DTKGLWidget::drawEverything()
{
   GLfloat white[] = {1.0f, 1.0f, 1.0f, 1.0f};
   glMaterialfv(GL_FRONT, GL_DIFFUSE, white);

   //drawAtoms();
   //drawLinks();
   drawCriticalPoints();
   drawGradientPaths();
}

void DTKGLWidget::addMolecule(QString fnam)
{
   QString wfnname=dtkglutils::getWFNFileNameFromCPX(fnam);

   gaussWaveFunc *wf_lcl=new gaussWaveFunc();
   if ( !(wf_lcl->readFromFile(wfnname.toStdString())) ) {
      QMessageBox::warning(this, tr("Error"),\
            (tr("Could not open the file ")+wfnname+tr("!")));
      delete wf_lcl;
      wf_lcl=NULL;
      return;
   }
   waveFunction.push_back(wf_lcl);

   DTKGLBondNetWork *bn_lcl=new DTKGLBondNetWork(this);
   if ( !(bn_lcl->readFromFile(wfnname)) ) {
      QMessageBox::warning(this, tr("Error"),\
            (tr("Could not open the file ")+wfnname+tr("!")));
      delete bn_lcl;
      bn_lcl=NULL;
      return;
   }
   if (cameraDistance<bn_lcl->getViewRadius()) {
      setCameraDistance(INITIAL_CAMERA_DISTANCE/bn_lcl->getViewRadius());
   }
   bondNW.push_back(bn_lcl);

   DTKGLCriticalPointNetWork *cp_lcl=new DTKGLCriticalPointNetWork(this);
   cp_lcl->setupRegularCPN(wf_lcl,bn_lcl);
   if ( !(cp_lcl->loadCPNStateFromFile(fnam)) ) {
      QMessageBox::warning(this, tr("Error"),\
            (tr("Could not open the file ")+fnam+tr("!")));
      delete cp_lcl;
      cp_lcl=NULL;
      return;
   }
   critPtNW.push_back(cp_lcl);
}

static void qNormalizeAngle(int &angle)
{
   while (angle < 0) angle += 360;
   while (angle >= 360) angle -= 360;
}

void DTKGLWidget::setXRotation(int angle)
{
   qNormalizeAngle(angle);
   if (angle != xRot) {
      xRot = angle;
      //emit xRotationChanged(angle);
      emit rotationChanged();
      update();
   }
}

void DTKGLWidget::setYRotation(int angle)
{
   qNormalizeAngle(angle);
   if (angle != yRot) {
      yRot = angle;
      //emit yRotationChanged(angle);
      emit rotationChanged();
      update();
   }
}

void DTKGLWidget::setZRotation(int angle)
{
   qNormalizeAngle(angle);
   if (angle != zRot) {
      zRot = angle;
      //emit zRotationChanged(angle);
      emit rotationChanged();
      update();
   }
}

void DTKGLWidget::setCameraDistance(double dist)
{
   if (dist!=cameraDistance) {
      cameraDistance=dist;
      emit zoomChanged();
      update();
   }
}

void DTKGLWidget::resetView()
{
  xRot=yRot=zRot=0.0f;
  emit rotationChanged();
  cameraDistance=INITIAL_CAMERA_DISTANCE;
  emit zoomChanged();
  update();
}

void DTKGLWidget::mousePressEvent(QMouseEvent *event)
{
   lastPos = event->pos();
}

void DTKGLWidget::mouseMoveEvent(QMouseEvent *event)
{
   int dx = event->x() - lastPos.x();
   int dy = event->y() - lastPos.y();

   if (event->buttons() & Qt::LeftButton) {
      setXRotation(xRot + dy/2 );
      setYRotation(yRot + dx/2 );
   } else if (event->buttons() & Qt::RightButton) {
      setXRotation(xRot + dy/2 );
      setZRotation(zRot + dx/2 );
   }

   lastPos = event->pos();
}

void DTKGLWidget::wheelEvent(QWheelEvent *event)
{
   int delta = event->delta();
   if (event->orientation() == Qt::Vertical) { if (delta < 0) {
         cameraDistance *= 0.8e0;
      } else if (delta > 0) {
         cameraDistance *= 1.25e0;
      }
      emit zoomChanged();
      update();
   }
   event->accept();
}

void DTKGLWidget::drawSingleSphere(float x, float y, float z, float radius,\
                                   float colr, float colg, float colb)
{
   glPushMatrix();
   glTranslatef(x,y,z);
   glColor3f(colr,colg,colb);
   GLUquadric *atomsph=gluNewQuadric();
   gluSphere(atomsph,radius,32,16);
   gluDeleteQuadric(atomsph);
   glPopMatrix();
}

void DTKGLWidget::drawSingleCylinder(QVector3D v0, float height, \
                                     float radius, float angle, QVector3D vrot, \
                                     QVector3D col)
{
   glPushMatrix();
   glTranslatef(v0[0], v0[1], v0[2]);
   glRotatef(angle,vrot[0],vrot[1],vrot[2]);
   glColor3f(col[0],col[1],col[2]);
   GLUquadric *cylinder=gluNewQuadric();
   gluCylinder(cylinder,radius,radius,height,24,1);
   gluDeleteQuadric(cylinder);
   glPopMatrix();
}



