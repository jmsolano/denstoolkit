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

    void setupRegularCPN(class gaussWaveFunc *uwf, DTKGLBondNetWork *ubn);
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

    bool iKnowACPs(void) {return cpn->iKnowACPs();}
    bool iKnowBCPs(void) {return cpn->iKnowBCPs();}
    bool iKnowRCPs(void) {return cpn->iKnowRCPs();}
    bool iKnowCCPs(void) {return cpn->iKnowCCPs();}
private:
    class critPtNetWork *cpn;
};

#endif // DTKGLCRITICALPOINTNETWORK_H