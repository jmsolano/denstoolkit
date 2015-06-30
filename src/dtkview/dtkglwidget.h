#ifndef DTKGLWIDGET_H
#define DTKGLWIDGET_H
#include <QOpenGLWidget>
#include <QMouseEvent>
#include <QMatrix4x4>
#include <QVector3D>

#define INITIAL_CAMERA_DISTANCE 2.5e0

/* This factor decreases the speed of rotation when a mouse drag is
 performed. The bigger the factor, the slower the change. */

class DTKGLWidget : public QOpenGLWidget
{
   Q_OBJECT
public:
   DTKGLWidget(QWidget *parent);
   void initializeGL();
   void paintGL();
   void resizeGL(int w,int h);
   void drawEverything();
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
   void setCameraDistance(double dist);
   void resetView(void);
signals:
   void rotationChanged(void);
   void xRotationChanged(int angle);
   void yRotationChanged(int angle);
   void zRotationChanged(int angle);
   void zoomChanged(void);
protected:
   void mousePressEvent(QMouseEvent *event);
   void mouseMoveEvent(QMouseEvent *event);
   void wheelEvent(QWheelEvent *event);
   void drawSingleSphere(float x,float y,float z,float radius,\
                         float colr,float colg,float colb);
   void drawSingleCylinder(QVector3D &v1, QVector3D &v2, \
                           float radius,
                           float colr, float colg, float colb);
   void getRotationVectorAndAngle(const QVector3D &v1,const QVector3D &v2,\
                                  QVector3D &vres,float &ares);
private:
   int xRot;
   int yRot;
   int zRot;
   double cameraDistance;
   QPoint lastPos;
   QMatrix4x4 pMatrix;
};

#endif // DTKGLWIDGET_H
