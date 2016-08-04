#ifndef _DEMAT1CRITPTNETWORKBP_H_
#define _DEMAT1CRITPTNETWORKBP_H_

/* ************************************************************************** */
class DeMat1CriticalPointNetworkBP {
/* ************************************************************************** */
public:
   DeMat1CriticalPointNetworkBP(class GaussWaveFunction &usrwf,\
         class critPtNetWork &usrcp);
   ~DeMat1CriticalPointNetworkBP();
   inline bool ImSetup(void) {return imsetup;}
   void SetupBondPath(int ata,int atb);
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   void init(void);
   void destroy(void);
   void SafetyChecks(void);
/* ************************************************************************** */
/* ************************************************************************** */
   class GaussWaveFunction *wf;
   class critPtNetWork *cp;
   int at1,at2;
/* ************************************************************************** */
private:
   bool imsetup;
   DeMat1CriticalPointNetworkBP(void) {} //Prohibited the use of default constructor.
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _DEMAT1CRITPTNETWORKBP_H_ */

