#ifndef DTKGLBONDNETWORK_H
#define DTKGLBONDNETWORK_H

#include <QObject>
#ifdef __APPLE__
#include <QOpenGLWidget>
#else
#include <QGLWidget>
#endif
#include <QVector3D>
#include "../common/bondnetwork.h"

class DTKGLBondNetWork : public QWidget
{
   Q_OBJECT
public:
   DTKGLBondNetWork(QWidget *parent = 0);
   ~DTKGLBondNetWork();
   bool readFromFile(QString filename);
   int numAtoms() {return atoms.size();}
   int numLinks() {return links.size();}
   QVector3D getAtomCoordinates(int idx) {return atoms[idx].r;}
   QVector3D getAtomColor(int idx) {return atoms[idx].color;}
   float getLinkAngle(int idx) {return links[idx].angle;}
   float getLinkHeight(int idx) {return links[idx].height;}
   QVector3D getLinkRotationVector(int idx) {return links[idx].rotVec;}
   QVector3D getLinkStart(int idx) {return links[idx].start;}
   QVector3D getLinkEnd(int idx) {return links[idx].end;}
   QVector3D getLinkColor(int idx) {return links[idx].color;}
   double getViewRadius(void) {if (bnw) {return ((bnw->rView)-(bnw->maxBondDist));}\
                               else {return 1.0e0;}}
   struct Link{
      float     angle;
      float     height;
      QVector3D rotVec;
      QVector3D start;
      QVector3D end;
      QVector3D color;
   };
   struct Atom{
      QVector3D r;
      QVector3D color;
   };
private:
   QVector<Link>   links;
   QVector<Atom>   atoms;
   bondNetWork *bnw;
};

#endif // DTKGLBONDNETWORK_H
