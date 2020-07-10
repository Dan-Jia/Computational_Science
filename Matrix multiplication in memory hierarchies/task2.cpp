// g++ task2.cpp
//./a.out 16 2 usw.
// ist b)ziemlich schneller als a), c)schneller als b)
#include <math.h>
#include <chrono>
#include <iostream>
#include <vector>

// ===================================
// Boost-Bibliothek für Alignment
#include <boost/align/aligned_allocator.hpp>
#include <boost/numeric/ublas/matrix.hpp>
using namespace boost::numeric::ublas;
// declare our 64-byte aligned allocator
template <typename T>
using aligned_alloc = boost::alignment::aligned_allocator<T, 64>;
// declare our aligned vector
template <typename T>
using aligned_vector = std::vector<T, aligned_alloc<T> >;

// a) ohne Kachelung
double skalProd;
void ohneKach(double* a, double* b, double* c, int n, int ld) {
  for (int i = 0; i < n; i++) {    // Zeilen
    for (int j = 0; j < n; j++) {  // Spalten
      skalProd = 0.0;
      for (int k = 0; k < n; k++) {
        skalProd += a[i * ld + k] * b[k * ld + j];
      }
      c[i * ld + j] += skalProd;
    }
  }
}

// b) Mit Kachelung
double *aBlock, *bBlock, *cBlock, *aRow, *bRow, *cRow;
void mitKach(double* a, double* b, double* c, int n, int ld, int m,
             int blockNum) {
  // loop Bloecke
  // Matrix c ist i*j; Matrix a ist i*k; Matrix b ist k*j
  for (int i = 0; i < blockNum; i++) {
    for (int k = 0; k < blockNum; k++) {
      for (int j = 0; j < blockNum; j++) {
        // Pointer auf erste Zeilenelemente der einzelnen Blöcke
        aBlock = a + (i * m) * ld + (k * m);
        bBlock = b + (k * m) * ld + (j * m);
        cBlock = c + (i * m) * ld + (j * m);

        // Enzelner Block zu berechnen
        for (int ii = 0; ii < m; ii++) {
          aRow = aBlock + ii * ld;
          cRow = cBlock + ii * ld;

          for (int kk = 0; kk < m; kk++) {
            bRow = bBlock + kk * ld;

            for (int jj = 0; jj < m; jj++) {
              cRow[jj] += aRow[kk] * bRow[jj];
            }
          }
        }
      }
    }
  }
}

// c) mit Blas-Routine DGEMM
extern "C" {
// C= alpha*A*B + beta*C
void dgemm_(char* TRANSA, char* TRANSB, const int* M, const int* N,
            const int* K, double* alpha, double* A, const int* LDA, double* B,
            const int* LDB, double* beta, double* C, const int* LDC);
}
char transXXX = 'N';
double alpha = 1.0;

// main Funktion
int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "error arguments" << std::endl;
    return -1;
  }

  int N = atoi(argv[1]);
  int m = atoi(argv[2]);

  int blockNum = ((int)sqrt(N) - 1) / m + 1;  // Anzahl der Bloecke
  int n = m * blockNum;  // Matrix Laenge mit dem ganzzahligen Fach von m

  // leading-dimension (64-byte-Alignment)
  // jede Zeile von Matrix gibt es ld Zahlen
  // jede cacheline gibt es 64/sizeof(double) Zahlen(double)
  // ((n - 1) / (64 / sizeof(double)) + 1) ist die Anzahl der cacheline
  // fuer eine Zeile der Matrix;
  int ld = ((n - 1) / (64 / sizeof(double)) + 1) * (64 / sizeof(double));

  // Matrixeintraege plus Padding. Standardmaessig mit Nullen intialisiert
  aligned_vector<double> a(ld * n);
  aligned_vector<double> b(n * ld);
  aligned_vector<double> c(n * n);

  int aAlign = (uintptr_t)a.data() % 64;
  int bAlign = (uintptr_t)b.data() % 64;
  int cAlign = (uintptr_t)c.data() % 64;
  int Padd = (ld * sizeof(double)) % 64;

  // Pruefen ob Matritzen aligned und gepadded sind (sollten alle aligned sein
  // aufgrund aligned_vector)
  if (aAlign != 0 || bAlign != 0 || cAlign != 0 || Padd != 0) {
    std::cout
        << std::endl
        << "---- ACHTUNG ---- Start-Adressen der Matrizen a, b, oder c sind "
           "nicht 64-byte aligned (Adresse ist kein Vielfaches von 64): "
        << aAlign << " " << bAlign << " " << cAlign << " " << Padd << std::endl
        << std::endl;
  }

  // Matrizen fuellen mit Zufallszahlen
  srand((unsigned)time(0));
  std::cout << "matrix a" << std::endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // kriegt die Zufallszahlen mit Groesse 1 bis 10
      a[i * ld + j] = 1 + (rand() % 10);
      std::cout << a[i * ld + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  std::cout << "matrix b" << std::endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // kriegt die Zufallszahlen mit Groesse 1 bis 10
      b[i * ld + j] = 1 + (rand() % 10);
      std::cout << b[i * ld + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  //=================================
  // Ergebnisse von a)
  //=================================
  int Iter = 1;
  using namespace std::chrono;
  auto start1 = system_clock::now();
  for (int i = 0; i < Iter; i++) ohneKach(a.data(), b.data(), c.data(), n, ld);
  auto end1 = system_clock::now();
  double time1 = duration<double>(end1 - start1).count();
  double gflops1 = 2 * pow(n, 3) * Iter * pow(10, -9) / time1;

  std::cout << "ohne Kachelung ist die Anzahl: " << gflops1 << "gFLOPS"
            << std::endl;

  // Ausgabe Matrix
  std::cout << "matrix c1" << std::endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      std::cout << c[i * ld + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  //===============================
  // Ergebnisse von b)
  //===============================
  fill(c.begin(), c.end(), 0);  // Matrix c zuruecksetzen
  auto start2 = system_clock::now();
  for (int i = 0; i < Iter; i++)
    mitKach(a.data(), b.data(), c.data(), n, ld, m, blockNum);
  auto end2 = system_clock::now();
  double time2 = duration<double>(end2 - start2).count();
  double gflops2 = 2 * pow(n, 3) * Iter * pow(10, -9) / time2;

  std::cout << "mit Kachelung ist die Anzahl: " << gflops2 << "gFLOPS"
            << std::endl;

  std::cout << "matrix c2" << std::endl;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      std::cout << c[i * ld + j] << " ";
    }
    std::cout << std::endl;
  }

  //===============================
  // Ergebnisse von c)
  //===============================
  fill(c.begin(), c.end(), 0);  // Matrix c zuruecksetzen
  auto start3 = system_clock::now();
  for (int i = 0; i < Iter; i++)
    dgemm_(&transXXX, &transXXX, &n, &n, &ld, &alpha, a.data(), &ld, b.data(),
           &ld, &alpha, c.data(), &ld);
  auto end3 = system_clock::now();
  double time3 = duration<double>(end3 - start3).count();
  double gflops3 = 2 * pow(n, 3) * Iter * pow(10, -9) / time3;

  std::cout << "mit der externe BLAS-Routine ist die Anzahl: " << gflops3
            << "gFLOPS" << std::endl;

  return 0;
}
