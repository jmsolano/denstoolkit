#ifndef _DEMAT1CRITPTNETWORKBP_H_
#define _DEMAT1CRITPTNETWORKBP_H_

/* ************************************************************************** */
class DeMat1CriticalPointNetworkBP {
   /* ************************************************************************** */
public:
   DeMat1CriticalPointNetworkBP(class GaussWaveFunction &usrwf,\
         class bondNetWork &usrbn);
   ~DeMat1CriticalPointNetworkBP();
   inline bool ImSetup(void) {return imsetup;}
   void ComputeCoreInteractionCPs2D(void);
   void ComputeCoreInteractionCPs6D(void);
   //void MapToUVCoordinatesM3x3(solreal (&e1)[3],solreal (&e2)[3],solreal (&huv)[2][2]);
   /* ************************************************************************** */
   class critPtNetWork *cpn;
   solreal *gCICP;
   solreal **hCICP;
   int nCICP;
   /* ************************************************************************** */
protected:
   /* ************************************************************************** */
   void init(void);
   void destroy(void);
   bool InitSafetyChecks(void);
   bool SetupCPN(void);
   void AllocAuxArrays(void);
   void ComputeSingleCICP2D(int idx);
   void ComputeSingleCICP6D(int idx);
   bool CPSafetyChecks(void);
   void GetTangentialVectors(const int bcpIdx,solreal (&e1)[3],solreal (&e2)[3]);
   int GetSignature(solreal (&v)[2]);
   int GetSignature(solreal (&v)[6]);
   /* ************************************************************************** */
   /* ************************************************************************** */
   class GaussWaveFunction *wf;
   class bondNetWork *bn;
   int at1,at2;
   /* ************************************************************************** */
private:
   bool imsetup;
   DeMat1CriticalPointNetworkBP(void) {} //Prohibited the use of default constructor.
   /* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _DEMAT1CRITPTNETWORKBP_H_ */

