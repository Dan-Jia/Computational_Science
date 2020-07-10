#include <math.h>
#include <time.h>
#include <x86intrin.h>
#include <iostream>

// direkte Berechnung der Kehrwerte
void dirKehr(float* x, float* res) {
  res[0] = 1 / x[0];
  res[1] = 1 / x[1];
  res[2] = 1 / x[2];
  res[3] = 1 / x[3];
}

// Berechnung der Kehrwerte mit SSE
void Ker_SSE(float* x, float* res) {
  __m128 s = _mm_setzero_ps();
  __m128 two = _mm_set_ps(2.0f, 2.0f, 2.0f, 2.0f);
  __m128 X = _mm_set_ps(x[3], x[2], x[1], x[0]);
  __m128 Y = _mm_rcp_ps(X);

  // Von der SSE-Nährungsfunktion bekommt man Startwerte Y der Newton Iteration
  // Newton Iterationsverfahren
  s = _mm_sub_ps(_mm_mul_ps(two, Y), _mm_mul_ps(X, _mm_mul_ps(Y, Y)));
  _mm_store_ps(&res[0], s);
}

// direkte Berechnung der Kehrwerte der Wurzel
void dieKehrWur(float* x, float* res) {
  res[0] = 1 / sqrt(x[0]);
  res[1] = 1 / sqrt(x[1]);
  res[2] = 1 / sqrt(x[2]);
  res[3] = 1 / sqrt(x[3]);
}

// Berechnung der Kehrwerte der Wurzel mit SSE
void KerWur_SSE(float* x, float* res) {
  __m128 s = _mm_setzero_ps();
  __m128 half = _mm_set_ps(0.5f, 0.5f, 0.5f, 0.5f);
  __m128 three_half = _mm_set_ps(1.5f, 1.5f, 1.5f, 1.5f);
  __m128 X = _mm_set_ps(x[3], x[2], x[1], x[0]);
  __m128 Y = _mm_rsqrt_ps(X);
  // Von der SSE-Nährungsfunktion bekommt man Startwerte Y der Newton Iteration
  // Newton Iterationsverfahren
  s = _mm_sub_ps(
      _mm_mul_ps(three_half, Y),
      _mm_mul_ps(half, _mm_mul_ps(X, _mm_mul_ps(Y, _mm_mul_ps(Y, Y)))));
  _mm_store_ps(res, s);
}

// main Funktion
int main(int argc, char* argv[]) {
  if (argc != 5) {
    std::cerr << "error arguments" << std::endl;
    return -1;
  }
  __m128* x;
  posix_memalign((void**)&x, 16, sizeof(float) * 4);

  ((float*)x)[0] = atoi(argv[1]);
  ((float*)x)[1] = atoi(argv[2]);
  ((float*)x)[2] = atoi(argv[3]);
  ((float*)x)[3] = atoi(argv[4]);

  float res[4];

  // Der Vergleich der Korrektheit und Genauigkeit
  clock_t start0 = clock();
  for (int i = 0; i < 100000; i++) dirKehr((float*)x, res);
  clock_t stop0 = (clock() - start0) * 1.0 / CLOCKS_PER_SEC * 1000;
  std::cout << "dirKehr für 100000 Wiederholungen: " << res[0] << " " << res[1]
            << " " << res[2] << " " << res[3] << ", " << stop0 << "ms"
            << std::endl;

  clock_t start1 = clock();
  for (int i = 0; i < 100000; i++) Ker_SSE((float*)x, res);
  clock_t stop1 = (clock() - start1) * 1.0 / CLOCKS_PER_SEC * 1000;
  std::cout << "Ker_SSE für 100000 Wiederholungen: " << res[0] << " " << res[1]
            << " " << res[2] << " " << res[3] << ", " << stop1 << "ms"
            << std::endl;

  clock_t start2 = clock();
  for (int i = 0; i < 100000; i++) dieKehrWur((float*)x, res);
  clock_t stop2 = (clock() - start2) * 1.0 / CLOCKS_PER_SEC * 1000;
  std::cout << "dieKehrWur für 100000 Wiederholungen: " << res[0] << " "
            << res[1] << " " << res[2] << " " << res[3] << ", " << stop2 << "ms"
            << std::endl;

  clock_t start3 = clock();
  for (int i = 0; i < 100000; i++) KerWur_SSE((float*)x, res);
  clock_t stop3 = (clock() - start3) * 1.0 / CLOCKS_PER_SEC * 1000;
  std::cout << "KerWur_SSE für 100000 Wiederholungen: " << res[0] << " "
            << res[1] << " " << res[2] << " " << res[3] << ", " << stop3 << "ms"
            << std::endl;

  return 0;
}
