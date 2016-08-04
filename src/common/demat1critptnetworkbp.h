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
   void SeekBondPath(int ata,int atb);
/* ************************************************************************** */
   class critPtNetWork *cpn;
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   void init(void);
   void destroy(void);
   bool SafetyChecks(void);
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

