// g++ task1.cpp -lpthread
// ./a.out
// die Ergebnisses:
// NUM_THREADS: 1, 2, 4, 8, 16, 32
// Zeit clock_gettime jeweils: 234.968 sec, 112.989 sec, 59.4477 sec,
// 54.2872 sec，49.2659 sec，51.904 sec

#include <pthread.h>
#include <stdlib.h>
#include <time.h>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>
#include <iostream>
#define maxIter 2000
using namespace boost::numeric::ublas;
// #define NUM_THREADS 1
// #define NUM_THREADS 2
// #define NUM_THREADS 4
// #define NUM_THREADS 8
// #define NUM_THREADS 16
#define NUM_THREADS 32

int N = 1024;  // Innere Gitterpunkte
matrix<double>*uNewPtr, *uOldPtr;
double a12, a21, a22, a23, a32;

// struct
struct messagetype {
  int n;    // Anteil pro Thread an Gesamtproblemgroesse n=N / NUM_THREADS
  int num;  // Die Nummer jedes Threads
};

// Jacobi-Verfahren
void* Jacobi(void* param) {
  struct messagetype* message = (messagetype*)param;

  // der Teil, der jedes Thread berechnen sollt.
  for (int i = (*message).num * (*message).n + 1;
       i < (*message).num * (*message).n + 1 + (*message).n; i++) {
    for (int j = 1; j < N + 1; j++) {
      (*uNewPtr)(i, j) =
          (1.0 - (*uOldPtr)(i - 1, j) * a21 - (*uOldPtr)(i + 1, j) * a23 -
           (*uOldPtr)(i, j - 1) * a12 - (*uOldPtr)(i, j + 1) * a32) /
          a22;
    }
  }
  return NULL;
}

// Die Vorbereitung des Jacobi-Verfahrens
void vorJacobiIter(matrix<double>& u, double vx, double vy) {
  const int N = u.size1();
  const double h = 1.0 / (N + 1.0);
  int r = 0;

  matrix<double> tempU(N + 2, N + 2);
  tempU = zero_matrix<double>(N + 2, N + 2);
  uNewPtr = &tempU;

  u = zero_matrix<double>(N + 2, N + 2);
  uOldPtr = &u;

  // linksseitige Differenzen durch den Fuenf-Punkt-Stern (A Matrix)
  a12 = -1.0 / h / h - vy / h;
  a21 = -1.0 / h / h - vx / h;
  a22 = 4.0 / h / h + vx / h + vy / h;
  a23 = -1.0 / h / h;
  a32 = -1.0 / h / h;

  // Iter (Vektor u)
  for (int iter = 0; iter < maxIter; iter++) {
    // Threads creating and join
    pthread_t tid, thread_id[NUM_THREADS];
    struct messagetype message[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; i++) {
      message[i].n = N / NUM_THREADS;
      message[i].num = i;

      r = pthread_create(&tid, NULL, Jacobi, (void*)&message[i]);
      thread_id[i] = tid;
    }
    for (int i = 0; i < NUM_THREADS; i++) {
      r = pthread_join(thread_id[i], NULL);
    }
    swap(*uOldPtr, *uNewPtr);  // uOldPtr und uNewPtr vertauschen
  }
}

// main Funktion
int main() {
  // die Parameter
  int v0 = 50;
  double vx = v0 / sqrt(5);
  double vy = 2.0 * v0 / sqrt(5);
  matrix<double> u(N + 2, N + 2);

  if (N % NUM_THREADS != 0) {
    std::cerr << "Die Vektoren können nicht gleichmäßig auf die Threads werden "
              << std::endl;
    return -1;
  }

  // Zeitmessung
  struct timespec start, stop;
  if (clock_gettime(CLOCK_REALTIME, &start) == -1) {
    perror("clock gettime");
    exit(EXIT_FAILURE);
  }

  // der Gebrauch des Iterationsverfahrens
  vorJacobiIter(u, vx, vy);

  if (clock_gettime(CLOCK_REALTIME, &stop) == -1) {
    perror("clock gettime");
    exit(EXIT_FAILURE);
  }
  double accum = ((stop.tv_sec - start.tv_sec) +
                  (stop.tv_nsec - start.tv_nsec) / 1000000000.0);

  int t = NUM_THREADS;
  std::cout << "NUM_THREADS: " << t << std::endl;
  std::cout << "Zeit clock_gettime: " << accum << " sec" << std::endl;

  // Ausgabe des Ergebnisses:
  for (int i = 0; i < N + 2; i++) {
    for (int j = 0; j < N + 2; j++) {
      std::cout << (*uOldPtr)(i, j) << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
