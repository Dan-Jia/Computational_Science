
#ifndef DISPLAY

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
}

void showVector(vector<float>* x, const int n) {
  for (int i = 0; i < n; i++) {
    printf("%1.3f \n", (*x)[i]);
  }  // i
  printf("\n");
}

#define DISPLAY
#endif
