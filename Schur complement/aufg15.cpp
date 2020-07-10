 #include <boost/chrono.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/thread.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

#include "display.cpp"
#include "schurmatrix.cpp"
#include "teilmatrix.cpp"

using namespace boost::numeric::ublas;

int main() {
  int n = 101;  // 41+2*Nullzeilen/-spalten
  TeilMatrix TM(n);
  SchurMatrix SM(&TM, n);

  vector<float>* f1 = new vector<float>(n * (n - 1) / 2, 1);
  vector<float>* f2 = new vector<float>(n * (n - 1) / 2, 1);
  vector<float>* fb = new vector<float>(n, 1);

  vector<float>* u1 = new vector<float>(n * (n - 1) / 2, 0);
  vector<float>* u2 = new vector<float>(n * (n - 1) / 2, 0);
  vector<float>* ub = new vector<float>(n, 0);

  vector<float>* z1 = new vector<float>(n * (n - 1) / 2, 0);
  vector<float>* z2 = new vector<float>(n * (n - 1) / 2, 0);
  vector<float>* z12 = new vector<float>(n, 0);
  vector<float>* z22 = new vector<float>(n, 0);

  vector<float>* t1 = new vector<float>(n * (n - 1) / 2, 0);
  vector<float>* t2 = new vector<float>(n * (n - 1) / 2, 0);

  CGSolverTM(&TM, 0, n * (n - 1) / 2, 0, n * (n - 1) / 2, z1, f1,
             n * (n - 1) / 2);  // z1=A11^-1*f1

  TM.mult(n * (n - 1) / 2, n, 0, n * (n - 1) / 2, z1,
          z12);  // z12=AB1*A11^-1*f1

  CGSolverTM(&TM, n * (n - 1) / 2 + n, n * (n - 1) / 2, n * (n - 1) / 2 + n,
             n * (n - 1) / 2, z2, f2,
             n * (n - 1) / 2);  // z2=A22^-1*f2

  TM.mult(n * (n - 1) / 2, n, n * (n - 1) / 2 + n, n * (n - 1) / 2, z2,
          z22);  // z22=A2B*A22^-1*f2

  (*fb) = (*fb) - (*z12) - (*z22);  // fb=fb(alt)+AB1*A11^-1*f1+A2B*A22^-1*f2

  CGSolverSM(&SM, ub, fb, n);  // ub=S^-1*fb

  TM.mult(0, n * (n - 1) / 2, n * (n - 1) / 2, n, ub,
          t1);            // t1=A1B*ub
  (*f1) = (*f1) - (*t1);  // f1=f1(alt)-A1B*ub

  TM.mult(n * (n - 1) / 2 + n, n * (n - 1) / 2, n * (n - 1) / 2, n, ub,
          t2);            // t2=AB2*ub
  (*f2) = (*f2) - (*t2);  // f2=f2(alt)-AB2*ub

  CGSolverTM(&TM, 0, n * (n - 1) / 2, 0, n * (n - 1) / 2, u1, f1,
             n * (n - 1) / 2);  // u1=A11^-1(f1-A1B*uB)

  CGSolverTM(&TM, n * (n - 1) / 2 + n, n * (n - 1) / 2, n * (n - 1) / 2 + n,
             n * (n - 1) / 2, u2, f2,
             n * (n - 1) / 2);  // u2=A22^-1(f2-AB2*uB)

  printf("\n");
  showMatrix(u1, (n - 1) / 2, n);
  showMatrix(ub, 1, n);
  showMatrix(u2, (n - 1) / 2, n);
  printf("\n");

  return 0;
}
