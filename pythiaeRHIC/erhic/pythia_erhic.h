#ifndef _PYTHIAERHIC_PYTHIA_ERHIC_H_
#define _PYTHIAERHIC_PYTHIA_ERHIC_H_

/*
 \file pythia_erhic.h
 Interfaces to entities in the eRHIC PYTHIA module.
 
 A note on name conversion between Fortran and C,
 true for gfortran at the time of writing.

 A subroutine or variable named x in Fortran
 will be named x_ (i.e. with a single trailing
 underscore) in C++.
 If the subroutine or variable is in a module named y,
 then the object x in Fortran will be named
 __y_MOD_x in C (i.e. TWO leading underscores, _MOD_
 between the module and object names, and no
 trailing underscore).
 Note that a common block declared in the module
 has the normal (i.e. single underscore) name, not
 the module convention.
 If in doubt, use nm to look for symbol names in
 the Fortran .o file.

 Also note that we put 'extern "C"' around all the
 variable and subroutine declarations to avoid
 C++ name mangling.
 
 \author Thomas Burton
 \date 10/9/12
 \copyright 2012 BNL
 */

// I'm usually averse to macros, but let's use one to generate
// the C++ equivalent variable names of Fortran variables.
// This should work for objects and subroutines:
//    FROM_FORTRAN_MODULE(pythia6, object)
//    FROM_FORTRAN_MODULE(pythia6, subroutine())
//#define FROM_FORTRAN_MODULE(module, variable) __##module##_MOD_##variable

extern "C" {
   // Initialisation subroutine.
   int __pythia6_MOD_initialise();
   int  __pythia6_MOD_generate();
   void __pythia6_MOD_finish();

   // Flag for printing ASCII output or not.
   // Set to 0 to suppress ASCII output, non-0 to generate it.
   extern int    __pythia6_MOD_printascii;
   extern int    __pythia6_MOD_genevent;
   extern double __pythia6_MOD_truex;
   extern double __pythia6_MOD_truew2;
   extern double __pythia6_MOD_truenu;
   extern double __pythia6_MOD_f1;
   extern double __pythia6_MOD_f2;
   extern double __pythia6_MOD_r;
   extern double __pythia6_MOD_sigma_rad;
   extern double __pythia6_MOD_sigradcor;
   extern double __pythia6_MOD_ebrems;
   extern double __pythia6_MOD_radgame;
   extern double __pythia6_MOD_radgamp[3];
} // extern "C"

#endif // _PYTHIAERHIC_PYTHIA_ERHIC_H_
