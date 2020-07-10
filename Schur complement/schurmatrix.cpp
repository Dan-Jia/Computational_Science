
#ifndef SCHURMATRIX

#include <boost/chrono.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/thread.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

#include "teilmatrix.cpp"

using namespace boost::numeric::ublas;

class SchurMatrix {
 public:
  SchurMatrix(TeilMatrix* TM, int n);
  int Sv(vector<float>* v, vector<float>* Sv);

 private:
  TeilMatrix* TM;
  int n;
  vector<float>* A1BuB;
  vector<float>* AB2uB;
  vector<float>* ub;
  vector<float>* s1;
  vector<float>* s2;
  vector<float>* t1;
  vector<float>* t2;
  vector<float>* t3;
};

SchurMatrix::SchurMatrix(TeilMatrix* TM, int n) {
  SchurMatrix::TM = TM;
  SchurMatrix::n = n;
  SchurMatrix::A1BuB = new vector<float>(n * (n - 1) / 2, 0);
  SchurMatrix::AB2uB = new vector<float>(n * (n - 1) / 2, 0);
  SchurMatrix::ub = new vector<float>(n, 1);
  SchurMatrix::s1 = new vector<float>(n * (n - 1) / 2, 0);
  SchurMatrix::s2 = new vector<float>(n * (n - 1) / 2, 0);
  SchurMatrix::t1 = new vector<float>(n, 0);
  SchurMatrix::t2 = new vector<float>(n, 0);
  SchurMatrix::t3 = new vector<float>(n, 0);
}

int SchurMatrix::Sv(vector<float>* v, vector<float>* Sv) {
  TM->mult(0, n * (n - 1) / 2, n * (n - 1) / 2, n, v, A1BuB);  // A1B*uB
  TM->mult(n * (n - 1) / 2 + n, n * (n - 1) / 2, n * (n - 1) / 2, n, v,
           AB2uB);  // AB2*uB
  CGSolverTM(TM, 0, n * (n - 1) / 2, 0, n * (n - 1) / 2, s1, A1BuB,
             n * (n - 1) / 2);  // A11^-1*A1B*uB
  CGSolverTM(TM, n * (n - 1) / 2 + n, n * (n - 1) / 2, n * (n - 1) / 2 + n,
             n * (n - 1) / 2, s2, AB2uB, n * (n + -1) / 2);  // A22^-1*AB2*uB
  TM->mult(n * (n - 1) / 2, n, 0, n * (n - 1) / 2, s1,
           t1);  // AB1*A11^-1*A1B*uB
  TM->mult(n * (n - 1) / 2, n, n * (n - 1) / 2 + n, n * (n - 1) / 2, s2,
           t2);                                             // A2B*A22^-1*AB2*uB
  TM->mult(n * (n - 1) / 2, n, n * (n - 1) / 2, n, v, t3);  // ABB*ub;
  (*Sv) = (*t3) - (*t1) - (*t2);
  return 0;
}

inline int CGSolverSM(SchurMatrix* SM, vector<float>* x, vector<float>* b,
                      const int n) {
  vector<float>* Ax = new vector<float>(n);
  vector<float>* r = new vector<float>(n);
  vector<float>* p = new vector<float>(n);
  vector<float>* Ap = new vector<float>(n, 0.0f);

  (*r) = (*b);                                  // r=b-A*x
  (*p) = (*r);                                  // p=r
  for (int iters = 0; iters < 1000; iters++) {  // Iteration
    SM->Sv(p, Ap);                              // A*p
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

#define SCHURMATRIX
#endif
