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
   int numAtoms() {if (bnw) {return bnw->nNuc;} else {return 0;}}
   QVector3D getAtomCoordinates(int idx) {return atoms[idx].r;}
   QVector3D getAtomColor(int idx) {return atoms[idx].color;}
   struct Link{
      float     angle;
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
