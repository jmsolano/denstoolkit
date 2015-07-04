#include "dtkglwidget.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <QOpenGLWidget>
#else
#include <GL/glut.h>
#include <QGLWidget>
#endif
#include <QtGui>
#include <QVector3D>
#include <QMatrix4x4>
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;
#include "dtkglutils.h"

#define DTKGL_DEFAULT_SPHERE_RADIUS 0.5f
#define DTKGL_DEFAULT_LINK_RADIUS 0.2f

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
   bondNW.clear();
   //zNear=0.01f;
   //zFar=100.0f;
}

DTKGLWidget::~DTKGLWidget()
{
   for(int i=0; i<bondNW.size(); ++i){
       delete bondNW[i];
   }
}

void DTKGLWidget::initializeGL()
{
   glClearColor(0.2,0.2,0.2,1.0);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_COLOR_MATERIAL);
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

void DTKGLWidget::drawEverything()
{
   GLfloat white[] = {1.0f, 1.0f, 1.0f, 1.0f};
   glMaterialfv(GL_FRONT, GL_DIFFUSE, white);

   drawAtoms();
   drawLinks();
}

void DTKGLWidget::addMolecule(QString fnam)
{
    DTKGLBondNetWork *bn_lcl=new DTKGLBondNetWork();
    bn_lcl->readFromFile(fnam);
    bondNW.push_back(bn_lcl);
    if (cameraDistance<bn_lcl->getViewRadius()) {
        setCameraDistance(INITIAL_CAMERA_DISTANCE/bn_lcl->getViewRadius());
    }
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



