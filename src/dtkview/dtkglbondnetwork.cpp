#include "dtkglbondnetwork.h"
#include <string>
#include "../common/bondnetwork.h"
#include "../common/atomcolschjmol.h"
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
   atoms.resize((bnw->nNuc));
   for (int i=0; i<atoms.size(); ++i) {
      for (int j=0; j<3; ++j) {atoms[i].r[j]=bnw->R[i][j];}
      atoms[i].color[0]=float(getAtomicRColorReal(bnw->atNum[i]));
      atoms[i].color[1]=float(getAtomicGColorReal(bnw->atNum[i]));
      atoms[i].color[2]=float(getAtomicBColorReal(bnw->atNum[i]));
   }
   links.resize(2*(bnw->nBonds));
   /*
   QVector3D va,vb;
   float alpha;
   int k=0,atni,atnk;
   solreal startpt[3],frak1;
   for (int i=0; i<nNuc; i++) {
      for (int j=0; j<MAXBONDINGATOMS; j++) {
         k=bNet[i][j];
         atni=atNum[i];
         atnk=atNum[k];
         //frak1=atomicRadius[atni]/(atomicRadius[atni]+atomicRadius[atnk]);
         frak1=getAtomicVDWRadius(atni)/(getAtomicVDWRadius(atni)+getAtomicVDWRadius(atnk));
         for (int l=0; l<3; l++) {
            startpt[l]=R[i][l]*(1.0e0-frak1)+R[k][l]*frak1;
         }
         if (k>0) {
            writePOVCylinder(pof,1,
                             R[i][0],R[i][1],R[i][2],
                  startpt[0],startpt[1],startpt[2],drawStickSize,
                  getAtomicRColorReal(atni),getAtomicGColorReal(atni),
                  getAtomicBColorReal(atni));
            writePOVCylinder(pof,1,
                             startpt[0],startpt[1],startpt[2],
                  R[k][0],R[k][1],R[k][2],drawStickSize,
                  getAtomicRColorReal(atnk),getAtomicGColorReal(atnk),
                  getAtomicBColorReal(atnk));
         }
      }
   }
   // */
   return true;
}
