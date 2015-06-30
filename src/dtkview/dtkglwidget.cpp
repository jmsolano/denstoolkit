#include "dtkglwidget.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <QVector3D>
#include <QMatrix4x4>
#include <cmath>
#include <iostream>
using std::cout;
using std::endl;

DTKGLWidget::DTKGLWidget(QWidget *parent)
   : QOpenGLWidget(parent)
{
   Q_UNUSED(parent);
   xRot=0;
   yRot=0;
   zRot=0;
   cameraDistance=INITIAL_CAMERA_DISTANCE;
   //zNear=0.01f;
   //zFar=100.0f;
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

void DTKGLWidget::drawEverything()
{
   GLfloat white[] = {1.0f, 1.0f, 1.0f, 1.0f};
   glMaterialfv(GL_FRONT, GL_DIFFUSE, white);

   QVector3D va(-1.0f,0.0f,0.0f);
   QVector3D vb(1.0f,0.0f,0.0f);

   float rad=0.5;
   drawSingleSphere(-1.0f,0.0f,0.0f, rad ,1.0f,0.6f,0.0f);
   drawSingleSphere(1.0f,0.0f,0.0f, rad ,1.0f,0.0f,0.1f);
   rad*=0.4f;
   drawSingleCylinder(va, vb, rad, 1.0f,0.6f,1.0f);
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

void DTKGLWidget::drawSingleCylinder(QVector3D &v1, QVector3D &v2, \
                                     float radius, \
                                     float colr, float colg, float colb)
{
   float height=v1.distanceToPoint(v2);
   float angle;
   QVector3D vrot;
   getRotationVectorAndAngle(v1,v2,vrot,angle);
   glPushMatrix();
   glTranslatef(v1[0], v1[1], v1[2]);
   glRotatef(angle,vrot[0],vrot[1],vrot[2]);
   glColor3f(colr,colg,colb);
   GLUquadric *cylinder=gluNewQuadric();
   gluCylinder(cylinder,radius,radius,height,24,1);
   gluDeleteQuadric(cylinder);
   glPopMatrix();
}

void DTKGLWidget::getRotationVectorAndAngle(const QVector3D &v1, const QVector3D &v2,\
                                            QVector3D &vres, float &ares)
{
    const float oeoPI = 180.0f/acos(-1.0f);
    static const QVector3D zAxis(0.0f, 0.0f, 1.0f);
    QVector3D diff=v2-v1;
    //float prod=QVector3D::dotProduct(diff,zAxis);
    vres=QVector3D::crossProduct(zAxis,diff);
    //float rad = acos(prod / ( v1.length() * v2.length()));
    float radians=acos(QVector3D::dotProduct(zAxis,diff)/(v1.length()*v2.length()));
    ares = radians * oeoPI;
}



