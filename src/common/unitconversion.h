#ifndef _UNITCONVERSION_H_
#define _UNITCONVERSION_H_

namespace unitconv {
   static constexpr double kCalPMole2Hartree=1.0e0/627.5e0;
   static constexpr double JPMole2Hartree=1.0e0/2.6255e+06;
   static constexpr double hartree2kJPerMole=2.6255e+03;
   static constexpr double hartree2JPerMole=2.6255e+06;
   static constexpr double hartree2kCalPerMole=6.27503e+02;
   static constexpr double hartree2cmm1=2.194746e+05;
   static constexpr double cmm12hartree=1.0e0/2.194746e+05;
   static constexpr double bohr2angstrom=0.529177249e0;
   static constexpr double angstrom2bohr=1.0e0/0.529177249e0;
}

#endif  /* _UNITCONVERSION_H_ */

