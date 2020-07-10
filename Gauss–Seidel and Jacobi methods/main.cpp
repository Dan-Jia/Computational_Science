#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <iostream>
using namespace boost::numeric::ublas;

// die Funktion von Jacobi-Verfahren
int jacobiIter(matrix<double>& u, const double vx, const double vy,
               const double a12, const double a21, const double a22,
               const double a23, const double a32) {
  const int n = u.size1();
  const double h = 1.0 / (n - 1.0);
  const int maxNumIter = 2000;
  matrix<double> tempU(n, n);
  tempU = zero_matrix<double>(n, n);
  u = zero_matrix<double>(n, n);
  matrix<double>* uNewPtr = &tempU;
  matrix<double>* uOldPtr = &u;
  matrix<double>* tempPtr;
  double maxFehler;
  int iter;
  int i;
  int j;

  for (iter = 0; iter < maxNumIter; iter++) {
    maxFehler = 0.0;
    for (i = 1; i < n - 1; i++) {
      for (j = 1; j < n - 1; j++) {
        (*uNewPtr)(i, j) = 1.0;  // =f
        (*uNewPtr)(i, j) -= (*uOldPtr)(i - 1, j) * a21;
        (*uNewPtr)(i, j) -= (*uOldPtr)(i + 1, j) * a23;
        (*uNewPtr)(i, j) -= (*uOldPtr)(i, j - 1) * a12;
        (*uNewPtr)(i, j) -= (*uOldPtr)(i, j + 1) * a32;
        (*uNewPtr)(i, j) = (*uNewPtr)(i, j) / a22;
        maxFehler =
            std::max(std::abs((*uNewPtr)(i, j) - (*uOldPtr)(i, j)), maxFehler);
      }                                   // i - Zeilen von u
    }                                     // j - Spalten von u
    if (maxFehler < 0.0000001) return 1;  // Genauigkeit erreicht
    tempPtr = uOldPtr;
    uOldPtr = uNewPtr;
    uNewPtr = tempPtr;
  }  // iter

  return 2;
}

// die Funktion von Gauss-Seidel-Verfahren
int gaussSeidelIter(matrix<double>& u, const double vx, const double vy,
                    const double a12, const double a21, const double a22,
                    const double a23, const double a32) {
  const int n = u.size1();
  const double h = 1.0 / (n - 1.0);
  const int maxNumIter = 1000000;
  u = zero_matrix<double>(n, n);
  double uIJOld;
  double maxFehler;
  int iter;
  int i;
  int j;

  for (iter = 0; iter < maxNumIter; iter++) {
    maxFehler = 0.0;
    for (i = 1; i < n - 1; i++) {
      for (j = 1; j < n - 1; j++) {
        uIJOld = u(i, j);
        u(i, j) = 1.0;  // =f
        u(i, j) -= u(i - 1, j) * a21;
        u(i, j) -= u(i + 1, j) * a23;
        u(i, j) -= u(i, j - 1) * a12;
        u(i, j) -= u(i, j + 1) * a32;
        u(i, j) = u(i, j) / a22;
        maxFehler = std::max(std::abs(u(i, j) - uIJOld), maxFehler);
      }                                   // i - Zeilen von u
    }                                     // j - Spalten von u
    if (maxFehler < 0.0000001) return 1;  // Genauigkeit erreicht
  }                                       // iter

  return 2;
}

// main Funktion
int main(int argc, char const* argv[]) {
  // die Parameter
  int n = 10;
  int v0 = 50;
  if (argc > 2) {
    n = atoi(argv[0]);
    v0 = atoi(argv[1]);
  }

  const double vx = v0 / sqrt(5);
  const double vy = 2.0 * v0 / sqrt(5);
  const double h = 1.0 / n;
  matrix<double> u(n + 1, n + 1);

  // zentrale Differenzen durch den Fuenf-Punkt-Ster
  const double a12 = -1.0 / (h * h) - vy / (2 * h);
  const double a21 = -1.0 / (h * h) - vx / (2 * h);
  const double a22 = 4.0 / (h * h);
  const double a23 = -1.0 / (h * h) + vx / (2 * h);
  const double a32 = -1.0 / (h * h) + vy / (2 * h);

  // linksseitige Differenzen durch den Fuenf-Punkt-Stern
  // const double a12 = -1.0 / h / h - vy / h;
  // const double a21 = -1.0 / h / h - vx / h;
  // const double a22 = 4.0 / h / h + vx / h + vy / h;
  // const double a23 = -1.0 / h / h;
  // const double a32 = -1.0 / h / h;

  // der Gebrauch des Iterationsverfahrens
  // int ret = jacobiIter(u, vx, vy, a12, a21, a22, a23, a32);
  int ret = gaussSeidelIter(u, vx, vy, a12, a21, a22, a23, a32);

  for (int i = 0; i < n + 1; i++) {
    for (int j = 0; j < n + 1; j++) {
      std::cout << u(i, j) << " ";
    }
    std::cout << std::endl;
  }
  /*
    // Ausgabe des Ergebnisses:
    for (int i = 0; i < n + 1; i++) {
      for (int j = 0; j < n + 1; j++) {
        std::cout << i * h << " " << j * h << " " << u(i, j) << std::endl;
      }
      std::cout << std::endl;
    }
  */
  return 1;
}
