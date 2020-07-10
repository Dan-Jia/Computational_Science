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

void showMatrix(vector<float>* x, const int n, const int m) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) printf("%1.3f \t", (*x)[i * m + j]);
    printf("\n");
  }  // i
  printf("\n");
}

void showResult(vector<float>* x1, vector<float>* x2, const int n, const int m,
                const int ovl) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m - ovl / 2; j++) printf("%1.3f \t", (*x1)[i * m + j]);
    for (int j = ovl / 2; j < m; j++) printf("%1.3f \t", (*x2)[i * m + j]);
    printf("\n");
  }  // i
  printf("\n");
}

inline void AMulP(vector<float>* p, vector<float>* Ap, const int n, const int m,
                  const float DivH2) {
  for (int i = 1; i < n - 1; i++) {
    for (int j = 1; j < m - 1; j++) {
      (*Ap)[i * m + j] =
          (4.0f * (*p)[i * m + j] - (*p)[i * m + j - 1] - (*p)[i * m + j + 1] -
           (*p)[(i - 1) * m + j] - (*p)[(i + 1) * m + j]) *
          DivH2;
    }  // j
  }    // i
}

inline int CGIter(vector<float>* x, vector<float>* p, vector<float>* r,
                  vector<float>* Ap, const int n, const int m,
                  const float DivH2) {
  AMulP(p, Ap, n, m, DivH2);               // A*p
  float pAp = inner_prod(*p, *Ap);         // p*A*p
  if (pAp < 0.00000001f) return -1;        // Abbruch wenn pAp zu klein
  float alpha = inner_prod(*p, *r) / pAp;  // alpha=p*r/(p*A*p)
  *x = *x + alpha * (*p);                  // x(neu)=x+alpha*p
  *r = *r - alpha * (*Ap);                 // r(neu)=r-alpha*A*p
  float beta = inner_prod(*Ap, *r) / pAp;  // beta=A*p*r/(A*p*p)
  *p = (*r) - beta * (*p);                 // p(neu)=r(neu)-beta*p
  return 0;
}

int main() {
  int n = 100;
  int m = 60;
  int ovl = (m - n / 2) * 2;
  int ohm2 = ovl - 1;
  int ohm1 = m - ovl;
  float H = 1.0f / ((float)n - 1.0f);
  float DivH2 = 1.0f / H / H;

  vector<float>* x1 = new vector<float>(n * m, 0.0f);
  vector<float>* x2 = new vector<float>(n * m, 0.0f);
  vector<float>* Ax1 = new vector<float>(n * m);
  vector<float>* Ax2 = new vector<float>(n * m);
  vector<float>* r1 = new vector<float>(n * m);
  vector<float>* r2 = new vector<float>(n * m);
  vector<float>* p1 = new vector<float>(n * m);
  vector<float>* p2 = new vector<float>(n * m);
  vector<float>* Ap1 = new vector<float>(n * m, 0.0f);
  vector<float>* Ap2 = new vector<float>(n * m, 0.0f);
  vector<float>* B1 = new vector<float>(n * m, 0.0f);
  vector<float>* B2 = new vector<float>(n * m, 0.0f);
  for (int i = 1; i < n - 1; i++)
    for (int j = 1; j < m - 1; j++) {
      (*B1)[i * m + j] = 1.0f;
      (*B2)[i * m + j] = 1.0f;
    }  // j

  for (int outerIter = 0; outerIter < 100; outerIter++) {
    for (int i = 1; i < n - 1; i++) {
      (*B1)[i * m + m - 2] = 1.0f + (*x2)[i * m + ohm2] * DivH2;
      (*B2)[i * m + 1] = 1.0f + (*x1)[i * m + ohm1] * DivH2;
    }  // i

    AMulP(x1, Ax1, n, m, DivH2);
    AMulP(x2, Ax2, n, m, DivH2);
    (*r1) = (*B1) - (*Ax1);  // r=b-A*x
    (*r2) = (*B2) - (*Ax2);  // r=b-A*x
    (*p1) = (*r1);           // p=r
    (*p2) = (*r2);           // p=r

    for (int iter = 0; iter < 1000; iter++) {
      if (CGIter(x1, p1, r1, Ap1, n, m, DivH2) == -1) break;
      if (CGIter(x2, p2, r2, Ap2, n, m, DivH2) == -1) break;
    }  // iter

  }  // outerIter

  // printf("ohm1: %i, ohm2: %i \n", ohm1, ohm2);
  // showMatrix(x1, n, m);
  // showMatrix(x2, n, m);
  // showResult(x1, x2, n, m, ovl);
  printf("u(0.5/0.5) = %1.5f \n", (*x1)[n / 2 * m + n / 2 - 1]);
  return 0;
}
