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

#ifndef DTKGLCRITICALPOINTNETWORK_H
#define DTKGLCRITICALPOINTNETWORK_H

#include <QWidget>
#include <QVector3D>
#include "dtkglbondnetwork.h"
#include "../common/critptnetwork.h"

class DTKGLCriticalPointNetWork : public QWidget
{
   Q_OBJECT
public:
    DTKGLCriticalPointNetWork(QWidget *parent = 0);
    ~DTKGLCriticalPointNetWork();

    void setupRegularCPN(class GaussWaveFunction *uwf, DTKGLBondNetWork *ubn);
    bool loadCPNStateFromFile(QString fname);

    int getNumACPs(void) {return cpn->nACP;}
    int getNumBCPs(void) {return cpn->nBCP;}
    int getNumRCPs(void) {return cpn->nRCP;}
    int getNumCCPs(void) {return cpn->nCCP;}
    int getNumPtsOfBGP(int bcpIdx) {return cpn->conBCP[bcpIdx][2];}
    int getBCPIdxInConnRCP(int rcpIdx,int bcpBox) {
       if (bcpBox<CPNW_MAXBCPSCONNECTEDTORCP) {
         return cpn->conRCP[rcpIdx][0][bcpBox];
       } else { return -1; }
    }
    int getNumPtsOfRGP(int rcpIdx,int bcpBox) {
       if (bcpBox<CPNW_MAXBCPSCONNECTEDTORCP) {
          return cpn->conRCP[rcpIdx][1][bcpBox];
       } else { return -1; }
    }
    int getRCPIdxInConnCCP(int ccpIdx,int rcpBox) {
       if (rcpBox<CPNW_MAXRCPSCONNECTEDTOCCP) {
         return cpn->conCCP[ccpIdx][0][rcpBox];
       } else { return -1; }
    }
    int getNumPtsOfCGP(int ccpIdx,int rcpBox) {
       if (rcpBox<CPNW_MAXRCPSCONNECTEDTOCCP) {
          return cpn->conCCP[ccpIdx][1][rcpBox];
       } else { return -1; }
    }
    QVector3D getACPCoordinates(int idx) {
       return QVector3D(cpn->RACP[idx][0],cpn->RACP[idx][1],cpn->RACP[idx][2]);
    }
    QVector3D getBCPCoordinates(int idx) {
       return QVector3D(cpn->RBCP[idx][0],cpn->RBCP[idx][1],cpn->RBCP[idx][2]);
    }
    QVector3D getRCPCoordinates(int idx) {
       return QVector3D(cpn->RRCP[idx][0],cpn->RRCP[idx][1],cpn->RRCP[idx][2]);
    }
    QVector3D getCCPCoordinates(int idx) {
       return QVector3D(cpn->RCCP[idx][0],cpn->RCCP[idx][1],cpn->RCCP[idx][2]);
    }

    QVector3D getBGPPointCoordinates(int bcpIdx,int ptIdx) {
       return QVector3D(cpn->RBGP[bcpIdx][ptIdx][0],
                        cpn->RBGP[bcpIdx][ptIdx][1],
                        cpn->RBGP[bcpIdx][ptIdx][2]);
    }
    QVector3D getRGPPointCoordinates(int rcpIdx,int bcpBox,int ptIdx){
       return QVector3D(cpn->RRGP[rcpIdx][bcpBox][ptIdx][0],
                        cpn->RRGP[rcpIdx][bcpBox][ptIdx][1],
                        cpn->RRGP[rcpIdx][bcpBox][ptIdx][2]);
    }
    QVector3D getCGPPointCoordinates(int ccpIdx,int rcpBox,int ptIdx){
       return QVector3D(cpn->RCGP[ccpIdx][rcpBox][ptIdx][0],
                        cpn->RCGP[ccpIdx][rcpBox][ptIdx][1],
                        cpn->RCGP[ccpIdx][rcpBox][ptIdx][2]);
    }
    QString getACPLabel(int idx) {return QString::fromStdString(cpn->lblACP[idx]);}
    QString getBCPLabel(int idx) {return QString::fromStdString(cpn->lblBCP[idx]);}
    QString getRCPLabel(int idx) {return QString::fromStdString(cpn->lblRCP[idx]);}
    QString getCCPLabel(int idx) {return QString::fromStdString(cpn->lblCCP[idx]);}

    bool iKnowACPs(void) {return cpn->IKnowACPs();}
    bool IKnowBCPs(void) {return cpn->IKnowBCPs();}
    bool iKnowRCPs(void) {return cpn->IKnowRCPs();}
    bool iKnowCCPs(void) {return cpn->IKnowCCPs();}
private:
    class CritPtNetWork *cpn;
};

#endif // DTKGLCRITICALPOINTNETWORK_H
