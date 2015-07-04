#include "dtkglbondnetwork.h"
#include <string>
#include "../common/bondnetwork.h"
#include "../common/atomcolschjmol.h"
#include "../common/atomradiicust.h"
#include <QMessageBox>
#include "dtkglutils.h"

DTKGLBondNetWork::DTKGLBondNetWork(QWidget *parent)
   : QWidget(parent)
{
   bnw=NULL;
   links.clear();
   atoms.clear();
}

DTKGLBondNetWork::~DTKGLBondNetWork()
{
   if(bnw!=NULL) {
       delete bnw;
       bnw=NULL;
   }
}

bool DTKGLBondNetWork::readFromFile(QString filename)
{
   std::string fnam=filename.toStdString();
   bnw=new bondNetWork();
   if ( !(bnw->readFromFile(fnam)) ) {
      QMessageBox::warning(this, tr("Error"),\
            tr("Could not open the file!"));
      delete bnw;
      bnw=NULL;
      return false;
   }
   bnw->setUpBNW();
   bnw->calcViewRadius();
   atoms.resize((bnw->nNuc));
   for (int i=0; i<atoms.size(); ++i) {
      for (int j=0; j<3; ++j) {atoms[i].r[j]=bnw->R[i][j];}
      atoms[i].color[0]=float(getAtomicRColorReal(bnw->atNum[i]));
      atoms[i].color[1]=float(getAtomicGColorReal(bnw->atNum[i]));
      atoms[i].color[2]=float(getAtomicBColorReal(bnw->atNum[i]));
   }
   links.resize(2*(bnw->nBonds));
   float alpha,frak1,dist;
   int k=0,atni,atnk,currIdx=0;
   QVector3D startpt,va,vb,rot;
   for (int i=0; i<bnw->nNuc; i++) {
      for (int j=0; j<MAXBONDINGATOMS; j++) {
         k=bnw->bNet[i][j];
         if (k>0) {
            atni=bnw->atNum[i];
            atnk=bnw->atNum[k];
            frak1=getAtomicVDWRadius(atni)/(getAtomicVDWRadius(atni)\
                                            +getAtomicVDWRadius(atnk));
            for (int l=0; l<3; l++) {
               va[l]=bnw->R[i][l];
               vb[l]=bnw->R[k][l];
               startpt[l]=va[l]*(1.0e0-frak1)+vb[l]*frak1;
            }
            dist=startpt.distanceToPoint(va);
            dtkglutils::getRotationVectorAndAngle(startpt,va,rot,alpha);
            links[currIdx].start=startpt;
            links[currIdx].end=va;
            links[currIdx].rotVec=rot;
            links[currIdx].angle=alpha;
            links[currIdx].height=dist;
            links[currIdx].color[0]=float(getAtomicRColorReal(atni));
            links[currIdx].color[1]=float(getAtomicGColorReal(atni));
            links[currIdx].color[2]=float(getAtomicBColorReal(atni));
            ++currIdx;
            dist=startpt.distanceToPoint(vb);
            dtkglutils::getRotationVectorAndAngle(startpt,vb,rot,alpha);
            links[currIdx].start=startpt;
            links[currIdx].end=vb;
            links[currIdx].rotVec=rot;
            links[currIdx].angle=alpha;
            links[currIdx].height=dist;
            links[currIdx].color[0]=float(getAtomicRColorReal(atnk));
            links[currIdx].color[1]=float(getAtomicGColorReal(atnk));
            links[currIdx].color[2]=float(getAtomicBColorReal(atnk));
            ++currIdx;
         }
      }
   }
   return true;
}
