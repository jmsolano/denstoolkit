#ifndef _DEMAT1CRITPTNETWORKBP_H_
#define _DEMAT1CRITPTNETWORKBP_H_

/* ************************************************************************** */
class DeMat1CriticalPointNetworkBP {
/* ************************************************************************** */
public:
   DeMat1CriticalPointNetworkBP(class GaussWaveFunction &usrwf,\
         class critPtNetWork &usrcp);
   ~DeMat1CriticalPointNetworkBP();
/* ************************************************************************** */
protected:
/* ************************************************************************** */
   void init();
   void destroy();
/* ************************************************************************** */
   DeMat1CriticalPointNetworkBP(void) {} //Prohibited the use of default constructor.
/* ************************************************************************** */
   class GaussWaveFunction *wf;
   class critPtNetWork *cp;
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _DEMAT1CRITPTNETWORKBP_H_ */

