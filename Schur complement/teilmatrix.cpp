#ifndef TEILMATRIX

#include <boost/chrono.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/thread.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

using namespace boost::numeric::ublas;

class TeilMatrix {
 public:
  TeilMatrix(int n);

  int mult(int kc, int nc, int kb, int nb, vector<float>* b, vector<float>* c);

 private:
  int n;
  float divH2;
};

TeilMatrix::TeilMatrix(int n) {
  TeilMatrix::n = n;
  TeilMatrix::divH2 = (n - 1.0f) * (n - 1.0f);
}

int TeilMatrix::mult(int kc, int nc, int kb, int nb, vector<float>* b,
                     vector<float>* c) {
  int iStart = kc;
  int iEnd = kc + nc - 1;
  int jStart = kb;
  int jEnd = kb + nb - 1;

  for (int i = iStart; i <= iEnd; i++) {
    float newVal = 0.0f;
    if ((jStart <= i) && (jEnd >= i)) {
      newVal += 4.0f * (*b)[i - jStart];
    }  // mult mit 4
    if ((jStart + 1 <= i) && (jEnd + 1 >= i) &&
        ((i - jStart - 1) % n != (n - 1))) {
      newVal += -1.0f * (*b)[i - jStart - 1];
    }  // mult mit -1
    if ((jStart <= i + 1) && (jEnd >= i + 1) && ((i - jStart + 1) % n != 0)) {
      newVal += -1.0f * (*b)[i - jStart + 1];
    }  // mult mit -1
    if ((jStart + n <= i) && (jEnd + n >= i) &&
        ((i - jStart - n) % n != (n - 1))) {
      newVal += -1.0f * (*b)[i - jStart - n];
    }  // mult mit -1
    if ((jStart <= i + n) && (jEnd >= i + n) && ((i - jStart + n) % n != 0)) {
      newVal += -1.0f * (*b)[i - jStart + n];
    }  // mult mit -1
    newVal *= divH2;
    (*c)[i - iStart] = newVal;
  }  // for i - Zeilen

  return 0;
}

inline int CGSolverTM(TeilMatrix* TM, int kc, int nc, int kb, int nb,
                      vector<float>* x, vector<float>* b, const int n) {
  vector<float>* Ax = new vector<float>(n);
  vector<float>* r = new vector<float>(n);
  vector<float>* p = new vector<float>(n);
  vector<float>* Ap = new vector<float>(n, 0.0f);

  (*r) = (*b);                                  // p=r=b-A*u
  (*p) = (*r);                                  // p=r=b-A*u
  for (int iters = 0; iters < 1000; iters++) {  // Iterationen
    TM->mult(kc, nc, kb, nb, p, Ap);            // A*p
    float pAp = inner_prod(*p, *Ap);            // p*A*p
    if (pAp < 0.000000001f) return -1;          // Abbruch wenn pAp zu klein
    float alpha = inner_prod(*p, *r) / pAp;     // alpha=p*r/(p*A*p)
    *x = *x + alpha * (*p);                     // x(neu)=x+alpha*p
    *r = *r - alpha * (*Ap);                    // r(neu)=r-alpha*A*p
    float beta = inner_prod(*Ap, *r) / pAp;     // beta=A*p*r/(A*p*p)
    *p = (*r) - beta * (*p);                    // p(neu)=r(neu)-beta*p
  }                                             // for iters
  return 0;
}

#define TEILMATRIX
#endif
