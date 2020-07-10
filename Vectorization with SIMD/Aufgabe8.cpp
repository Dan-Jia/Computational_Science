#include <math.h>
#include <time.h>
#include <x86intrin.h>
#include <iostream>

// skale Berechnung der Werte mit der exp Funktion
void skaleExp(float* x, float* res) {
  res[0] = exp(x[0]);
  res[1] = exp(x[1]);
  res[2] = exp(x[2]);
  res[3] = exp(x[3]);
}

// Berechnung mit SSE-Vektorperationen(mit fester Zahl)
void festZahl_SSE(float* x, int n, float* res) {
  __m128 X = _mm_set_ps(x[3], x[2], x[1], x[0]);
  __m128 one = _mm_set_ps(1.0f, 1.0f, 1.0f, 1.0f);
  __m128 s = _mm_set_ps(1.0f, 1.0f, 1.0f, 1.0f);

  for (int i = n; i > 0; i--) {
    __m128 I_recip = _mm_set_ps(1.0f / (float)i, 1.0f / (float)i,
                                1.0f / (float)i, 1.0f / (float)i);
    s = _mm_add_ps(_mm_mul_ps(s, _mm_mul_ps(X, I_recip)), one);
    // s= s*(X/i)+1
  }
  _mm_store_ps(&res[0], s);
}

// main Funktion
int main(int argc, char* argv[]) {
  if (argc != 6) {
    std::cerr << "error arguments" << std::endl;
    return -1;
  }

  int n = atoi(argv[1]);
  __m128* x;
  posix_memalign((void**)&x, 16, sizeof(float) * 4);
  ((float*)x)[0] = atoi(argv[2]);
  ((float*)x)[1] = atoi(argv[3]);
  ((float*)x)[2] = atoi(argv[4]);
  ((float*)x)[3] = atoi(argv[5]);

  float res[4];

  // Der Vergleich der Korrektheit und Genauigkeit
  clock_t start0 = clock();
  for (int i = 0; i < 100000; i++) skaleExp((float*)x, res);
  clock_t stop0 = (clock() - start0) * 1.0 / CLOCKS_PER_SEC * 1000;
  std::cout << "skaleExp für 100000 Wiederholungen: " << res[0] << " " << res[1]
            << " " << res[2] << " " << res[3] << ", " << stop0 << "ms"
            << std::endl;

  clock_t start1 = clock();
  for (int i = 0; i < 100000; i++) festZahl_SSE((float*)x, n, res);
  clock_t stop1 = (clock() - start1) * 1.0 / CLOCKS_PER_SEC * 1000;
  std::cout << "festZahl_SSE für 100000 Wiederholungen: " << res[0] << " "
            << res[1] << " " << res[2] << " " << res[3] << ", " << stop1 << "ms"
            << std::endl;

  return 0;
}
