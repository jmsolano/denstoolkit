/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.6.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
                             Santiago Alberto Flores Roman
          Copyright (c) 2013-2022, Juan Manuel Solano-Altamirano
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
#include <QPainter>
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
#define DTKGL_DEFAULT_TRANSPARENCY 0.4f
#define DTKGL_DISPLACEMENT_CP_LABEL_X 2
#define DTKGL_DISPLACEMENT_CP_LABEL_Y 2

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
   drawAts=setTransp=true;
   drawAtLbls=drawBnds=false;
   drawBGPs=drawRGPs=drawCGPs=true;
   drawCPLbls=false;
   cameraDistance=INITIAL_CAMERA_DISTANCE;
   waveFunction.clear();
   bondNW.clear();
   critPtNW.clear();
   //zNear=0.01f;
   //zFar=100.0f;
}
DTKGLWidget::~DTKGLWidget() {
   clearWFsBNsCPXs();
}
void DTKGLWidget::clearWFsBNsCPXs() {
   for (int i=0; i<waveFunction.size(); ++i) {
      delete waveFunction[i];
      waveFunction[i]=NULL;
   }
   waveFunction.clear();
   for (int i=0; i<bondNW.size(); ++i){
      delete bondNW[i];
      bondNW[i]=NULL;
   }
   bondNW.clear();
   for (int i=0; i<critPtNW.size(); ++i) {
      delete critPtNW[i];
      critPtNW[i]=NULL;
   }
   critPtNW.clear();
}
void DTKGLWidget::initializeGL() {
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
   //glClearColor((48.f/255.f),(62.f/255.f),(115.f/255.f),1.0f); // MediumBlue01
   glClearColor((57.f/255.f),(105.f/255.f),(196.f/255.f),1.0f); // Ocean
   glShadeModel (GL_SMOOTH);

   glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);

   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_COLOR_MATERIAL);
   glEnable(GL_NORMALIZE); //Needed for keeping light uniform
   glEnable(GL_BLEND); //Enable blending.
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //Set blending function.
}
void DTKGLWidget::paintGL() {

   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glLoadIdentity();
   glTranslatef(0.0, 0.0, -10.0);
   glRotatef(float(xRot), 1.0, 0.0, 0.0);
   glRotatef(float(yRot), 0.0, 1.0, 0.0);
   glRotatef(float(zRot), 0.0, 0.0, 1.0);
   glScalef(cameraDistance,cameraDistance,cameraDistance);

   drawEverything();
}
void DTKGLWidget::resizeGL(int w, int h) {
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
void DTKGLWidget::drawAtoms() {
   DTKGLBondNetWork *bn;
   QVector3D pos,col;
   float rad = DTKGL_DEFAULT_SPHERE_RADIUS;
   if (setTransp) {
      for (int i=0; i<bondNW.size(); ++i) {
         bn=bondNW[i];
         for (int j=0; j<(bn->numAtoms()); ++j) {
            pos=bn->getAtomCoordinates(j);
            col=bn->getAtomColor(j);
            drawSingleTransparentSphere(pos,rad,col,DTKGL_DEFAULT_TRANSPARENCY);
         }
      }
   } else {
      for (int i=0; i<bondNW.size(); ++i) {
         bn=bondNW[i];
         for (int j=0; j<(bn->numAtoms()); ++j) {
            pos=bn->getAtomCoordinates(j);
            col=bn->getAtomColor(j);
            drawSingleSphere(
            pos[0],pos[1],pos[2],
            rad, col[0], col[1], col[2]
            );
         }
      }
   }
}
void DTKGLWidget::drawLinks() {
    DTKGLBondNetWork *bn;
    float rad = DTKGL_DEFAULT_LINK_RADIUS;
    if (setTransp) {
       for (int bIdx=0; bIdx<bondNW.size(); ++bIdx) {
          bn=bondNW[bIdx];
          for (int i=0; i<bn->numLinks(); ++i) {
             drawSingleTransparentCylinder(\
                      bn->getLinkStart(i),bn->getLinkHeight(i),rad,\
                      bn->getLinkAngle(i),bn->getLinkRotationVector(i),\
                      bn->getLinkColor(i),DTKGL_DEFAULT_TRANSPARENCY);
          }
       }
    } else {
       for (int bIdx=0; bIdx<bondNW.size(); ++bIdx) {
          bn=bondNW[bIdx];
          for (int i=0; i<bn->numLinks(); ++i) {
             drawSingleCylinder(bn->getLinkStart(i),bn->getLinkHeight(i),rad,\
                      bn->getLinkAngle(i),bn->getLinkRotationVector(i),\
                      bn->getLinkColor(i));
          }
       }
    }
}
void DTKGLWidget::drawAtomLabels() {
   DTKGLBondNetWork *bn;
   QVector3D pos;
   QString lbl;
   for (int bIdx=0; bIdx<bondNW.size(); ++bIdx) {
      bn=bondNW[bIdx];
      for (int i=0; i<(bn->numAtoms()); ++i) {
         lbl=bn->getAtomLabel(i);
         pos=bn->getAtomCoordinates(i);
         drawText(pos,lbl);
      }
   }
}
void DTKGLWidget::drawAttractorCriticalPoints() {
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
void DTKGLWidget::drawBondCriticalPoints() {
   DTKGLCriticalPointNetWork *cp;
   const static QVector3D colBCP(0.0f,0.6f,1.0f);
   QVector3D pos;
   float rad=DTKGL_DEFAULT_CRIT_PT_RADIUS;
   for (int cpIdx=0; cpIdx<critPtNW.size(); ++cpIdx) {
      cp=critPtNW[cpIdx];
      if (!(cp->IKnowBCPs())) {return;}
      for (int i=0; i<cp->getNumBCPs(); ++i) {
         pos=cp->getBCPCoordinates(i);
         drawSingleSphere(pos,rad,colBCP);
      }
   }
}
void DTKGLWidget::drawRingCriticalPoints() {
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
void DTKGLWidget::drawCageCriticalPoints() {
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
void DTKGLWidget::drawCriticalPoints() {
   drawAttractorCriticalPoints();
   drawBondCriticalPoints();
   drawRingCriticalPoints();
   drawCageCriticalPoints();
}
void DTKGLWidget::drawCPLabels() {
   DTKGLCriticalPointNetWork *cpn;
   QVector3D pos;
   QString lbl;
   int nn;
   for (int cIdx=0; cIdx<critPtNW.size(); ++cIdx) {
      cpn=critPtNW[cIdx];
      nn=cpn->getNumACPs();
      for (int i=0; i<nn; ++i) {
         lbl=tr("A")+QString::number(i+1);
         pos=cpn->getACPCoordinates(i);
         drawText(pos,lbl,DTKGL_DISPLACEMENT_CP_LABEL_X,DTKGL_DISPLACEMENT_CP_LABEL_Y);
      }
      nn=cpn->getNumBCPs();
      for (int i=0; i<nn; ++i) {
         lbl=tr("B")+QString::number(i+1);
         pos=cpn->getBCPCoordinates(i);
         drawText(pos,lbl,DTKGL_DISPLACEMENT_CP_LABEL_X,DTKGL_DISPLACEMENT_CP_LABEL_Y);
      }
      nn=cpn->getNumRCPs();
      for (int i=0; i<nn; ++i) {
         lbl=tr("R")+QString::number(i+1);
         pos=cpn->getRCPCoordinates(i);
         drawText(pos,lbl,DTKGL_DISPLACEMENT_CP_LABEL_X,DTKGL_DISPLACEMENT_CP_LABEL_Y);
      }
      nn=cpn->getNumCCPs();
      for (int i=0; i<nn; ++i) {
         lbl=tr("C")+QString::number(i+1);
         pos=cpn->getCCPCoordinates(i);
         drawText(pos,lbl,DTKGL_DISPLACEMENT_CP_LABEL_X,DTKGL_DISPLACEMENT_CP_LABEL_Y);
      }
   }
}
void DTKGLWidget::drawBondGradientPaths() {
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
void DTKGLWidget::drawRingGradientPaths() {
   DTKGLCriticalPointNetWork *cp;
   const static QVector3D colRGP(0.0f,0.8f,0.0f);
   QVector3D pos;
   int npts,k;
   float rad=DTKGL_DEFAULT_GRAD_PATH_RADIUS;
   for (int cpIdx=0; cpIdx<critPtNW.size(); ++cpIdx) {
      cp=critPtNW[cpIdx];
      for (int i=0; i<cp->getNumRCPs(); ++i) {
         k=0;
         while (cp->getBCPIdxInConnRCP(i,k)>=0) {
            npts=cp->getNumPtsOfRGP(i,k);
            for (int l=0; l<npts; ++l) {
               pos=cp->getRGPPointCoordinates(i,k,l);
               drawSingleSphere(pos,rad,colRGP);
            }
            ++k;
         }
      }
   }
}
void DTKGLWidget::drawCageGradientPaths() {
   DTKGLCriticalPointNetWork *cp;
   const static QVector3D colCGP(1.0f,0.5f,0.0f);
   QVector3D pos;
   int npts,k;
   float rad=DTKGL_DEFAULT_GRAD_PATH_RADIUS;
   for (int cpIdx=0; cpIdx<critPtNW.size(); ++cpIdx) {
      cp=critPtNW[cpIdx];
      for (int i=0; i<cp->getNumCCPs(); ++i) {
         k=0;
         while (cp->getRCPIdxInConnCCP(i,k)>=0) {
            npts=cp->getNumPtsOfCGP(i,k);
            for (int l=0; l<npts; ++l) {
               pos=cp->getCGPPointCoordinates(i,k,l);
               drawSingleSphere(pos,rad,colCGP);
            }
            ++k;
         }
      }
   }
}
void DTKGLWidget::drawGradientPaths() {
   if (drawBGPs) {drawBondGradientPaths();}
   if (drawRGPs) {drawRingGradientPaths();}
   if (drawCGPs) {drawCageGradientPaths();}
}
void DTKGLWidget::drawEverything() {
   GLfloat white[] = {1.0f, 1.0f, 1.0f, 1.0f};
   glMaterialfv(GL_FRONT, GL_DIFFUSE, white);

   drawCriticalPoints();
   drawGradientPaths();
   if (drawAts) { drawAtoms(); }
   if (drawAtLbls) { drawAtomLabels(); }
   if (drawCPLbls) { drawCPLabels(); }
   if (drawBnds) { drawLinks(); }
}
void DTKGLWidget::addMolecule(QString fnam) {
   QString wfnname=dtkglutils::getWFNFileNameFromCPX(fnam);

   GaussWaveFunction *wf_lcl=new GaussWaveFunction();
   if ( !(wf_lcl->ReadFromFile(wfnname.toStdString())) ) {
      QMessageBox::warning(this, tr("Error"),\
            (tr("Could not open the file ")+wfnname+tr("!")\
             +tr("\nPlease remember that cpx and its associated"\
                 " wfn(wfx) file must be in the same directory.")));
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
static void qNormalizeAngle(int &angle) {
   while (angle < 0) angle += 360;
   while (angle >= 360) angle -= 360;
}
void DTKGLWidget::setXRotation(int angle) {
   qNormalizeAngle(angle);
   if (angle != xRot) {
      xRot = angle;
      //emit xRotationChanged(angle);
      emit rotationChanged();
      update();
   }
}
void DTKGLWidget::setYRotation(int angle) {
   qNormalizeAngle(angle);
   if (angle != yRot) {
      yRot = angle;
      //emit yRotationChanged(angle);
      emit rotationChanged();
      update();
   }
}
void DTKGLWidget::setZRotation(int angle) {
   qNormalizeAngle(angle);
   if (angle != zRot) {
      zRot = angle;
      //emit zRotationChanged(angle);
      emit rotationChanged();
      update();
   }
}
void DTKGLWidget::setViewAtoms(bool val) {
   drawAts=val;
   update();
}
void DTKGLWidget::setDrawAtomLabels(bool dal) {
   drawAtLbls=dal;
   emit drawAtLblsChanged();
   update();
}
void DTKGLWidget::setViewRegularBonds(bool vrb) {
   drawBnds=vrb;
   update();
}
void DTKGLWidget::setViewBondGradientPaths(bool dbgp) {
   drawBGPs=dbgp;
   update();
}
void DTKGLWidget::setViewRingGradientPaths(bool drgp) {
   drawRGPs=drgp;
   update();
}
void DTKGLWidget::setViewCageGradientPaths(bool dcgp) {
   drawCGPs=dcgp;
   update();
}
void DTKGLWidget::setDrawCPLabels(bool dal) {
   drawCPLbls=dal;
   update();
}
void DTKGLWidget::setTransparentAtomsAndLinks(bool val) {
   setTransp=val;
   update();
}
void DTKGLWidget::setCameraDistance(double dist) {
   if (dist!=cameraDistance) {
      cameraDistance=dist;
      emit zoomChanged();
      update();
   }
}
void DTKGLWidget::resetView() {
  xRot=yRot=zRot=0.0f;
  emit rotationChanged();
  cameraDistance=INITIAL_CAMERA_DISTANCE;
  emit zoomChanged();
  update();
}
void DTKGLWidget::mousePressEvent(QMouseEvent *event) {
   lastPos = event->pos();
}
void DTKGLWidget::mouseMoveEvent(QMouseEvent *event) {
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
void DTKGLWidget::wheelEvent(QWheelEvent *event) {
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
                                   float colr, float colg, float colb) {
   glPushMatrix();
   glTranslatef(x,y,z);
   glColor3f(colr,colg,colb);
   GLUquadric *atomsph=gluNewQuadric();
   gluSphere(atomsph,radius,16,8);
   gluDeleteQuadric(atomsph);
   glPopMatrix();
}
void DTKGLWidget::drawSingleTransparentSphere(QVector3D r, float rad, QVector3D c, float t) {
   glPushMatrix();
   glTranslatef(r[0],r[1],r[2]);
   glColor4f(c[0],c[1],c[2],t);
   GLUquadric *atomsph=gluNewQuadric();
   gluSphere(atomsph,rad,32,16);
   gluDeleteQuadric(atomsph);
   glPopMatrix();
}
void DTKGLWidget::drawSingleCylinder(QVector3D v0, float height, \
                                     float radius, float angle, QVector3D vrot, \
                                     QVector3D col) {
   glPushMatrix();
   glTranslatef(v0[0], v0[1], v0[2]);
   glRotatef(angle,vrot[0],vrot[1],vrot[2]);
   glColor3f(col[0],col[1],col[2]);
   GLUquadric *cylinder=gluNewQuadric();
   gluCylinder(cylinder,radius,radius,height,24,1);
   gluDeleteQuadric(cylinder);
   glPopMatrix();
}
void DTKGLWidget::drawSingleTransparentCylinder(QVector3D v0, float height, float radius, float angle, QVector3D vrot, QVector3D col, float t) {
   glPushMatrix();
   glTranslatef(v0[0], v0[1], v0[2]);
   glRotatef(angle,vrot[0],vrot[1],vrot[2]);
   glColor4f(col[0],col[1],col[2],t);
   GLUquadric *cylinder=gluNewQuadric();
   gluCylinder(cylinder,radius,radius,height,24,1);
   gluDeleteQuadric(cylinder);
   glPopMatrix();
}
void DTKGLWidget::drawText(QVector3D r, QString lbl,int dx,int dy) {
   GLdouble modelview[16],projection[16];
   GLint viewport[4];
   GLdouble wX,wY,wZ;
   glGetDoublev(GL_MODELVIEW_MATRIX,modelview);
   glGetDoublev(GL_PROJECTION_MATRIX, projection);
   glGetIntegerv(GL_VIEWPORT,viewport);
   GLdouble x,y,z;
   x=GLdouble(r[0]);
   y=GLdouble(r[1]);
   z=GLdouble(r[2]);
   gluProject(x,y,z,modelview,projection,viewport,&wX,&wY,&wZ);
   wX+=dx; wY+=dy;

   glMatrixMode(GL_PROJECTION);
   glPushMatrix();
   glLoadIdentity();
   gluOrtho2D(0.0, float(this->width()), 0.0, float(this->height()));

   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();

   glDisable(GL_LIGHTING);
   glColor3f(1.0, 1.0, 1.0); // Green
   glRasterPos2i(GLint(wX),GLint(wY));
   string s=lbl.toStdString();
   void * font = GLUT_BITMAP_HELVETICA_12;
   for (string::iterator i = s.begin(); i != s.end(); ++i)
   {
       glutBitmapCharacter(font, *i);
   }
   glEnable(GL_LIGHTING);
   glMatrixMode(GL_PROJECTION);
   glPopMatrix();

   glMatrixMode(GL_MODELVIEW);
   glPopMatrix();
}


