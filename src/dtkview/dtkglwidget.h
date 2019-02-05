/*
 *                    This source code is part of
 *
 *                  D  E  N  S  T  O  O  L  K  I  T
 *
 *                         VERSION: 1.3.1
 *
 *             Contributors: Juan Manuel Solano-Altamirano
 *                           Julio Manuel Hernandez-Perez
 *        Copyright (c) 2013-2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 * If you want to redistribute modifications of the suite, please
 * consider to include your modifications in our official release.
 * We will be pleased to consider the inclusion of your code
 * within the official distribution. Please keep in mind that
 * scientific software is very special, and version control is 
 * crucial for tracing bugs. If in despite of this you distribute
 * your modified version, please do not call it DensToolKit.
 *
 * If you find DensToolKit useful, we humbly ask that you cite
 * the paper(s) on the package --- you can find them on the top
 * README file.
 *
 */

/* Several aspects of this class were inspired by the program
 * VizMol, whose author is Abraham Max Santos Ramos.
 */

#ifndef DTKGLWIDGET_H
#define DTKGLWIDGET_H
#ifdef __APPLE__
#include <QOpenGLWidget>
#else
#include <QGLWidget>
#endif
#include <QMouseEvent>
#include <QMatrix4x4>
#include <QtGui>
#include <QVector3D>
#include "../common/gausswavefunction.h"
#include "dtkglbondnetwork.h"
#include "dtkglcriticalpointnetwork.h"

#define INITIAL_CAMERA_DISTANCE 2.5e0

/* This factor decreases the speed of rotation when a mouse drag is
 performed. The bigger the factor, the slower the change. */
#ifdef __APPLE__
class DTKGLWidget : public QOpenGLWidget
#else
class DTKGLWidget : public QGLWidget
#endif
{
   Q_OBJECT
public:
   DTKGLWidget(QWidget *parent);
   ~DTKGLWidget();
   void clearWFsBNsCPXs(void);
   void initializeGL();
   void paintGL();
   void resizeGL(int w,int h);
   void drawAtoms();
   void drawLinks();
   void drawAtomLabels();
   void drawAttractorCriticalPoints();
   void drawBondCriticalPoints();
   void drawRingCriticalPoints();
   void drawCageCriticalPoints();
   void drawCriticalPoints();
   void drawCPLabels();
   void drawBondGradientPaths();
   void drawRingGradientPaths();
   void drawCageGradientPaths();
   void drawGradientPaths();
   void drawEverything();
   void addMolecule(QString fnam);
   int getXRot(void) {return (xRot);}
   int getYRot(void) {return (yRot);}
   int getZRot(void) {return (zRot);}
   double getCameraDistance(void) {return cameraDistance;}
   int getCurrentZoom(void) {return int(100*cameraDistance/INITIAL_CAMERA_DISTANCE);}
public slots:
   // slots for xyz-rotation slider
   void setXRotation(int angle);
   void setYRotation(int angle);
   void setZRotation(int angle);
   void setViewAtoms(bool val);
   void setDrawAtomLabels(bool dal);
   void setViewRegularBonds(bool vrb);
   void setViewBondGradientPaths(bool dbgp);
   void setViewRingGradientPaths(bool drgp);
   void setViewCageGradientPaths(bool dcgp);
   void setDrawCPLabels(bool dal);
   void setTransparentAtomsAndLinks(bool val);
   void setCameraDistance(double dist);
   void resetView(void);
signals:
   void rotationChanged(void);
   void xRotationChanged(int angle);
   void yRotationChanged(int angle);
   void zRotationChanged(int angle);
   void zoomChanged(void);
   void drawAtLblsChanged(void);
protected:
   void mousePressEvent(QMouseEvent *event);
   void mouseMoveEvent(QMouseEvent *event);
   void wheelEvent(QWheelEvent *event);
   void drawSingleSphere(float x,float y,float z,float radius,\
                         float colr,float colg,float colb);
   void drawSingleSphere(QVector3D r,float rad,QVector3D c) {
      return drawSingleSphere(r[0],r[1],r[2],rad,c[0],c[1],c[2]);
   }
   void drawSingleTransparentSphere(QVector3D r,float rad,QVector3D c,float t);
   void drawSingleCylinder(QVector3D v0, float height, \
                           float radius, float angle,QVector3D vrot,
                           QVector3D col);
   void drawSingleTransparentCylinder(QVector3D v0, float height, \
                           float radius, float angle,QVector3D vrot,
                           QVector3D col,float t);
   void drawText(QVector3D r, QString lbl, int dx=0, int dy=0);
private:
   int xRot;
   int yRot;
   int zRot;
   bool drawAts;
   bool drawAtLbls;
   bool drawBnds;
   bool drawBGPs;
   bool drawRGPs;
   bool drawCGPs;
   bool drawCPLbls;
   bool setTransp;
   double cameraDistance;
   QPoint lastPos;
   QMatrix4x4 pMatrix;
   QVector<GaussWaveFunction*> waveFunction;
   QVector<DTKGLBondNetWork*> bondNW;
   QVector<DTKGLCriticalPointNetWork*> critPtNW;
};

#endif // DTKGLWIDGET_H
