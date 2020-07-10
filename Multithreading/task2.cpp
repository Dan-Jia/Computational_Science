// g++ task2.cpp -lpthread
// ./a.out
// die Ergebnisses:
// NUM_THREADS: 1, 2, 4, 8, 16, 32
// Zeit clock_gettime jeweils: 347.172 sec，301.584 sec，268.406 sec，
// 288.105 sec, 288.656 sec, 287.657 sec

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
// int Laenge;
float skalar_rr, skalar_ar;
vector<float>*r, *a, *u;
pthread_mutex_t mtx;

// struct
struct messagetype {
  int n;    // Anteil pro Thread an Gesamtproblemgroesse
  int num;  // Die Nummer jedes Threads
};

// Skalarprodukte mit pthreads
void* skalarPro(void* param) {
  struct messagetype* message = (messagetype*)param;
  float teilSkalar_rr = 0.0f;
  float teilSkalar_ar = 0.0f;

  for (int i = (*message).num * (*message).n;
       i < (*message).num * (*message).n + (*message).n; i++) {
    teilSkalar_ar += (*a)[i] * (*r)[i];
    teilSkalar_rr += (*r)[i] * (*r)[i];
  }

  pthread_mutex_lock(&mtx);
  skalar_ar += teilSkalar_ar;
  skalar_rr += teilSkalar_rr;
  pthread_mutex_unlock(&mtx);

  return NULL;
}

// Die Vorbereitung der Skalarprodukte
void vorSkalarPro() {
  // Threads creating and join
  int m = 0;
  pthread_t tid, thread_id[NUM_THREADS];
  struct messagetype message[NUM_THREADS];
  pthread_mutex_init(&mtx, NULL);

  for (int i = 0; i < NUM_THREADS; i++) {
    // Anteil pro Thread an Gesamtproblemgroesse
    message[i].n = (N + 2) * (N + 2) / NUM_THREADS;
    message[i].num = i;  // Der Nummer jedes Threads

    m = pthread_create(&tid, NULL, skalarPro, (void*)&message[i]);
    thread_id[i] = tid;
  }

  for (int i = 0; i < NUM_THREADS; i++) {
    m = pthread_join(thread_id[i], NULL);
  }
}

// Gradientenverfahren
void Gradient(float vx, float vy) {
  float h = 1.0 / (N + 1.0);
  int Laenge = N + 2;
  int Laenge_vek = (N + 2) * (N + 2);

  a = new vector<float>(Laenge_vek, 0);  // a = A*r
  u = new vector<float>(Laenge_vek, 0);  // u = u + alpha * r
  r = new vector<float>(Laenge_vek, 0);  // r = r - alpha * a

  // linksseitige Differenzen durch den Fuenf-Punkt-Stern (A Matrix)
  float a12 = -1.0 / h / h - vy / h;
  float a21 = -1.0 / h / h - vx / h;
  float a22 = 4.0 / h / h + vx / h + vy / h;
  float a23 = -1.0 / h / h;
  float a32 = -1.0 / h / h;

  // Startwerte fue r0 (wenn u0=0), r0 = b-A*u0 = b = 1
  for (int i = Laenge; i < Laenge_vek - Laenge; i++)
    if (i % Laenge != 0 && i % Laenge != Laenge - 1) {
      (*r)[i] = 1;
    }

  // Iter (Vektor a = A*r)
  for (int iter = 0; iter < maxIter; iter++) {
    for (int i = 1; i < Laenge - 1; i++) {
      for (int j = i * Laenge + 1; j < (i + 1) * Laenge - 1; j++) {
        (*a)[j] = a21 * (*r)[j - 1] + a12 * (*r)[j - Laenge] +
                  a23 * (*r)[j + 1] + a32 * (*r)[j + Laenge] + a22 * (*r)[j];
      }
    }

    // u(k+1) und r(k+1) berechnen
    skalar_rr = 0.0f;
    skalar_ar = 0.0f;
    vorSkalarPro();
    float alpha = skalar_rr / skalar_ar;

    for (int i = 0; i < Laenge_vek; i++) {
      (*u)[i] = (*u)[i] + alpha * (*r)[i];  // u = u + alpha * r
      (*r)[i] = (*r)[i] - alpha * (*a)[i];  // r = r - alpha * a
    }
  }
}

// main Funktion
int main() {
  // die Parameter
  int v0 = 50;
  float vx = v0 / sqrt(5);
  float vy = 2.0 * v0 / sqrt(5);

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
  Gradient(vx, vy);

  if (clock_gettime(CLOCK_REALTIME, &stop) == -1) {
    perror("clock gettime");
    exit(EXIT_FAILURE);
  }
  float accum = ((stop.tv_sec - start.tv_sec) +
                 (stop.tv_nsec - start.tv_nsec) / 1000000000.0);

  int t = NUM_THREADS;
  std::cout << "NUM_THREADS: " << t << std::endl;
  std::cout << "Zeit clock_gettime: " << accum << " sec" << std::endl;

  // Ausgabe des Ergebnisses:
  for (int i = 0; i < N + 2; i++) {
    for (int j = 0; j < N + 2; j++) {
      std::cout << (*u)(j * (N + 2) + i) << " ";
    }
    std::cout << std::endl;
  }

  free(a);
  free(r);
  free(u);

  return 0;
}
