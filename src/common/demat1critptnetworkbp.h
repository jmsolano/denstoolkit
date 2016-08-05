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
   void ComputeCoreInteractionCPs(void);
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
   void ComputeSingleCICP(int idx);
   bool CPSafetyChecks(void);
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

