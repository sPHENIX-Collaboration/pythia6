/**
 \file pythia_commons.cpp
 
 \author Thomas Burton 
 \date 10/9/12
 \copyright 2012 BNL.
 */

#include "pythia_commons.h"

int msti(int i) {
   return pypars_.msti[i - 1];
}

double pari(int i) {
   return pypars_.pari[i - 1];
}

double vint(int i) {
   return pyint1_.vint[i - 1];
}

int k(int i, int j) {
   return pyjets_.k[j - 1][i - 1];
}

double p(int i, int j) {
   return pyjets_.p[j - 1][i - 1];
}

double v(int i, int j) {
   return pyjets_.v[j - 1][i - 1];
}

// Note that the first array index runs [0:500]
// unlike a default Fortran array, so this index is NOT decremented
// when accessed in C++. Array indices are reversed as usual so
// ngen(0, 3) in Fortran becomes ngen[2, 0] in C++.
int ngen(int i, int j) {
   return pyint5_.ngen[j - 1][i];
}