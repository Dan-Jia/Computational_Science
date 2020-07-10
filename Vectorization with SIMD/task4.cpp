#include <math.h>
#include <time.h>
#include <x86intrin.h>
#include <iostream>

// Jacobi-Iteration ohne SIMD Vektoren
void jacobiIter(const int n, const float vx, const float vy, const float h,
                float* u, const float a12, const float a21, const float a22,
                const float a23, const float a32) {
  const int maxNumIter = 2000;

  for (long int Iter = 0; Iter < maxNumIter; Iter++) {
    for (long int i = 1; i < n; i++) {
      for (long int j = 1; j < n; j++) {
        u[i * n + j] =
            (1.0 - u[i * n + (j - 1)] * a21 - u[(i - 1) * n + j] * a12 -
             u[(i + 1) * n + j] * a32 - u[i * n + (j + 1)] * a23) /
            a22;
      }
    }
  }
}

// Jacobi-Iteration mit SIMD Vektoren
void jacobi_SSE(const int n, const float vx, const float vy, const float h,
                __m128* u_SSE, const float a12, const float a21,
                const float a22, const float a23, const float a32) {
  const int maxNumIter = 2000;
  const float tol = 0.0000000001f;
  float* uOld;
  uOld = (float*)u_SSE;
  __m128 A12 = _mm_set1_ps(a12);
  __m128 A32 = _mm_set1_ps(a32);
  __m128 A21 = _mm_set1_ps(a21);
  __m128 A23 = _mm_set1_ps(a23);
  __m128 A22 = _mm_set1_ps(a22);
  __m128 one = _mm_set1_ps(1.0f);
  __m128 zero = _mm_setzero_ps();
  __m128 Tol = _mm_set1_ps(tol);

  // Vektorisierung
  for (long int Iter = 0; Iter < maxNumIter; Iter++) {
    for (long int i = 1; i < n; i++) {
      for (long int j = 1; j < n / 4 + 1; j++) {
        __m128 leftU = _mm_set_ps(
            uOld[i * (n + 8) + j * 4 + 2], uOld[i * (n + 8) + j * 4 + 1],
            uOld[i * (n + 8) + j * 4], uOld[i * (n + 8) + j * 4 - 1]);
        __m128 rightU = _mm_set_ps(
            uOld[i * (n + 8) + j * 4 + 4], uOld[i * (n + 8) + j * 4 + 3],
            uOld[i * (n + 8) + j * 4 + 2], uOld[i * (n + 8) + j * 4 + 1]);

        __m128 leftProd = _mm_mul_ps((leftU), (A21));
        __m128 rightProd = _mm_mul_ps((rightU), (A23));
        __m128 upperProd = _mm_mul_ps(u_SSE[(i - 1) * (n / 4 + 2) + j], (A12));
        __m128 lowerProd = _mm_mul_ps(u_SSE[(i + 1) * (n / 4 + 2) + j], (A32));
        // Iterationsschritt:
        u_SSE[i * (n / 4 + 2) + j] = _mm_div_ps(
            _mm_sub_ps(
                _mm_sub_ps(_mm_sub_ps(_mm_sub_ps(one, upperProd), lowerProd),
                           leftProd),
                rightProd),
            A22);
      }
    }
  }
}

// main Funktion
int main(int argc, char const* argv[]) {
  if (argc != 2) {
    std::cerr << "error arguments" << std::endl;
    return -1;
  }
  // die Parameter
  int v0 = atoi(argv[1]);
  const int n = 1024;
  const float vx = v0 / sqrt(5.0f);
  const float vy = v0 * 2.0f / sqrt(5.0f);
  const float h = 1.0f / (n + 1.0f);

  float* u = new float[(n + 2) * (n + 2)];
  __m128* u_SSE;
  posix_memalign((void**)&u_SSE, 16, sizeof(float) * (n + 2) * (n + 8));

  // linksseitige Differenzen durch den Fuenf-Punkt-Stern
  const float a12 = -1.0f / (h * h) - vy / h;
  const float a21 = -1.0f / (h * h) - vx / h;
  const float a22 = 1.0f / (4.0f / (h * h) + vx / h + vy / h);
  const float a23 = -1 / (h * h);
  const float a32 = -1 / (h * h);

  // Der Vergleich der Korrektheit
  clock_t start0 = clock();
  jacobiIter(n, vx, vy, h, u, a12, a21, a22, a23, a32);
  clock_t stop0 = (clock() - start0) * 1.0 / CLOCKS_PER_SEC * 1000;
  std::cout << "Berechnungsdauer von jacobiIter: " << stop0 << "ms"
            << std::endl;

  clock_t start1 = clock();
  jacobi_SSE(n, vx, vy, h, u_SSE, a12, a21, a22, a23, a32);
  clock_t stop1 = (clock() - start1) * 1.0 / CLOCKS_PER_SEC * 1000;
  std::cout << "Berechnungsdauer von jacobi_SSE: " << stop1 << "ms"
            << std::endl;

  return 0;
}
