#ifndef _INTEGRATEOVERBONDPATH_H_
#define _INTEGRATEOVERBONDPATH_H_

/* ************************************************************************** */
class IntegrateOverBondPath {
/* ************************************************************************** */
public:
   IntegrateOverBondPath(class critPtNetWork &ucpn);
/* ************************************************************************** */
protected:
   class critPtNetWork *cp;
   /* Functions  */
   void init(void);
/* ************************************************************************** */
};
/* ************************************************************************** */


#endif  /* _INTEGRATEOVERBONDPATH_H_ */

