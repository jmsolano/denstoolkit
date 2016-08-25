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
   solreal **eivalCICP2D; /*!< Contains the eigenvalues of 2D CICPs (projected gamma)  */
   solreal **eivalCICP6D; /*!< Contains the eigenvalues of 6D CICPs (projected gamma)  */
   int *sigCICP2D; /*!< Contains the information of 2D signatures (projected gamma)  */
   int *sigCICP6D; /*!< Contains the information of 6D signatures.  */
   int nCICP;
   /* ************************************************************************** */
protected:
   /* ************************************************************************** */
   void init(void);
   void destroy(void);
   bool InitSafetyChecks(void);
   bool SetupCPN(void);
   bool AllocAuxArrays(void);
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

