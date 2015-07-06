#include "dtkglcriticalpointnetwork.h"
#include <QDebug>
#include "../common/wavefunctionclass.h"
#include "../common/bondnetwork.h"
#include "../common/critptnetwork.h"

DTKGLCriticalPointNetWork::DTKGLCriticalPointNetWork(QWidget *parent)
      : QWidget(parent)
{
   cpn=NULL;
}

DTKGLCriticalPointNetWork::~DTKGLCriticalPointNetWork()
{
  if (cpn!=NULL) {
     delete cpn;
     cpn=NULL;
  }
}

void DTKGLCriticalPointNetWork::setupRegularCPN(gaussWaveFunc *uwf, DTKGLBondNetWork *ubn)
{
   cpn=new critPtNetWork(*uwf,*(ubn->bnw));

}

bool DTKGLCriticalPointNetWork::loadCPNStateFromFile(QString fname)
{
   if (cpn==NULL) {
      qDebug() << "Error: first load the critical point network!";
      return false;
   }
   return cpn->readFromFile(fname.toStdString());
}

