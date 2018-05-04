#ifndef _PYTHIAERHIC_PYTHIA_COMMONS_H_
#define _PYTHIAERHIC_PYTHIA_COMMONS_H_

/**
 \file pythia_commons.h
 Interfaces to PYTHIA common blocks
 Remember that Fortran arrays are indexed the opposite
 way round to C arrays, and array bounds are different by default:
 e.g. K(1,100) in Fortran becomes k[99][0] in C++.
 
 \author Thomas Burton
 \date 10/9/12
 \copyright 2012 BNL
 */

extern "C" {
   // COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
   extern struct {
      int n;
      int npad;
      int k[5][4000];
      double p[5][4000];
      double v[5][4000];
   } pyjets_;
   
   // COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
   extern struct {
      int mstp[200];
      double parp[200];
      int msti[200];
      double pari[200];
   } pypars_;
   
   // COMMON/PYINT1/MINT(400),VINT(400)
   extern struct {
      int mint[400];
      double vint[400];
   } pyint1_;
   
   // COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
   extern struct {
      int ngenpd;
      int ngen[3][501];
      double xsec[3][501];
   } pyint5_;
} // extern "C"

// =====================================================================
// Wrapper functions around PYTHIA common block variables.
// Takes care of different Fortran and C array indexing
// so they can be called with the PYTHIA syntax e.g. k(1, 2).
// =====================================================================

int msti(int);

double pari(int);

double vint(int);

int k(int, int);

double p(int, int);

double v(int, int);

int ngen(int, int);

#endif // _PYTHIAERHIC_PYTHIA_COMMONS_H_
