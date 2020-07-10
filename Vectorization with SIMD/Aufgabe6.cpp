#include <time.h>
#include <x86intrin.h>
#include <iostream>

// das Skalarprodukt ohne SIMD-Verwendung
float normal_skalPro(const float* a, const float* b, int n) {
  float res = 0;

  for (int i = 0; i < n; ++i) {
    res += a[i] * b[i];
  }
  return res;
}

// das Skalarprodukt mit SIMD-Verwendung
float simd_skalPro(const float* a, const float* b, int n) {
  __m128 A, B;
  __m128 s = _mm_setzero_ps();
  // Diese Varialbe s wird gegeben und initialisiert, um die Ergebnisse von A+B
  // zu speichern.
  float res = 0.0, temp[4];

  for (int i = 0; i + 4 < n; i += 4) {
    A = _mm_load_ps(a + i);
    B = _mm_load_ps(b + i);

    s = _mm_add_ps(s, _mm_mul_ps(A, B));
    // Jede schleife macht die Summen von je 4 Komponenten von A und
    // B(a0*b0,a1*b1,a2*b2,a3*b3),
    // Am Ende bekommt man 4 Summen,danach werden die 4 Summen in s gespeichert.
  }
  _mm_store_ps(&temp[0], s);
  // laden 4 Summen s (s0,s1,s2,s3) im Hauptspeicher
  res = temp[0] + temp[1] + temp[2] + temp[3];

  return res;
}

// main Funktion
int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "error arguments" << std::endl;
    return -1;
  }

  int n = atoi(argv[1]);
  float* a;
  float* b;
  posix_memalign((void**)&a, 16, sizeof(float) * n);
  posix_memalign((void**)&b, 16, sizeof(float) * n);

  for (int i = 0; i < n; ++i) {
    a[i] = i + 1;
    b[i] = n - i;
  }

  float res0, res1;
  clock_t start0 = clock();
  for (int i = 0; i < 100000; i++) res0 = normal_skalPro(a, b, n);
  clock_t stop0 = (clock() - start0) * 1.0 / CLOCKS_PER_SEC * 1000;
  std::cout << "normal_skalPro für 100000 Wiederholungen: " << res0 << ", "
            << stop0 << "ms" << std::endl;

  clock_t start1 = clock();
  for (int i = 0; i < 100000; i++) res1 = simd_skalPro(a, b, n);
  clock_t stop1 = (clock() - start1) * 1.0 / CLOCKS_PER_SEC * 1000;
  std::cout << "simd_skalPro für 100000 Wiederholungen: " << res1 << ", "
            << stop1 << "ms" << std::endl;

  free(a);
  free(b);
  return 0;
}
