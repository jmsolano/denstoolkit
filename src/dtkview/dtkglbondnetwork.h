/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
                           VERSION: 1.5.0
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
                             Luis Alfredo Nunez-Meneses
          Copyright (c) 2013-2021, Juan Manuel Solano-Altamirano
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

#ifndef DTKGLBONDNETWORK_H
#define DTKGLBONDNETWORK_H

#include <QObject>
#ifdef __APPLE__
#include <QOpenGLWidget>
#else
#include <QGLWidget>
#endif
#include <QVector3D>
#include <QString>
#include "../common/bondnetwork.h"

class DTKGLBondNetWork : public QWidget
{
   Q_OBJECT
   friend class DTKGLCriticalPointNetWork;
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
   QString getAtomLabel(int idx) {return QString::fromStdString(bnw->atLbl[idx]);}
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
protected:
   BondNetWork *bnw;
private:
   QVector<Link>   links;
   QVector<Atom>   atoms;
};

#endif // DTKGLBONDNETWORK_H
