/*
                      This source code is part of
  
                    D  E  N  S  T  O  O  L  K  I  T
  
               Contributors: Juan Manuel Solano-Altamirano
                             Julio Manuel Hernandez-Perez
          Copyright (c) 2013-2020, Juan Manuel Solano-Altamirano
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

#ifndef _DEMAT1CRITPTNETWORKBP_H_
#define _DEMAT1CRITPTNETWORKBP_H_

/* ************************************************************************** */
class DeMat1CriticalPointNetworkBP {
   /* ************************************************************************** */
public:
   DeMat1CriticalPointNetworkBP(class GaussWaveFunction &usrwf,\
         class BondNetWork &usrbn);
   ~DeMat1CriticalPointNetworkBP();
   inline bool ImSetup(void) {return imsetup;}
   void ComputeCoreInteractionCPs2D(void);
   void ComputeCoreInteractionCPs6D(void);
   bool DifferentSignaturesCICPvsNN(void);
   //void MapToUVCoordinatesM3x3(double (&e1)[3],double (&e2)[3],double (&huv)[2][2]);
   class CritPtNetWork *cpn;
   double **eivalCICP2D; /*!< Contains the eigenvalues of 2D CICPs (projected gamma, ACP-ACP connection)  */
   double **eivalCICP6D; /*!< Contains the eigenvalues of 6D CICPs (ACP-ACP connection) */
   double **eivalNN6D; /*!< Contains the eigenvalues of 6D Nuc-Nuc correlations  */
   int *sigCICP2D; /*!< Contains the information of 2D cicp signatures (projected gamma)  */
   int *sigCICP6D; /*!< Contains the information of 6D cicp signatures.  */
   int *sigNN6D; /*!< Contains the information of 6D nuc-nuc signatures.  */
   int nCICP;
   /* ************************************************************************** */
protected:
   /* ************************************************************************** */
   void Init(void);
   void Destroy(void);
   bool InitSafetyChecks(void);
   bool SetupCPN(void);
   bool AllocAuxArrays(void);
   void ComputeSingleCICP2D(int idx);
   void ComputeSingleCICP6D(int idx);
   bool CPSafetyChecks(void);
   void GetTangentialVectors(const int bcpIdx,double (&e1)[3],double (&e2)[3]);
   int GetSignature(double (&v)[2]);
   int GetSignature(double (&v)[6]);
   /* ************************************************************************** */
   class GaussWaveFunction *wf;
   class BondNetWork *bn;
   int at1,at2;
   /* ************************************************************************** */
protected:
   bool imsetup;
   DeMat1CriticalPointNetworkBP(void) {} //Prohibited the use of default constructor.
   void AssignHessian6D(double (&hh)[3][3],double (&hph)[3][3],double (&hp)[3][3],\
         double (&hess)[6][6]);
   /* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _DEMAT1CRITPTNETWORKBP_H_ */

