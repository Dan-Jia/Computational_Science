#include <x86intrin.h>
#include <algorithm>
#include <cmath>
#include <iostream>

void showMatrix(__m128* M128, int n, int m) {  // zum Auslese nach Matlab z.B.
  float* M = (float*)M128;
  std::cout << std::fixed;
  std::cout << std::endl;
  for (int i = 0; i < n; i++) {        // Zeilen
    for (int j = 3; j < m - 3; j++) {  // Spalten
      std::cout << M[i * m + j] << " ";
    }  // j - Spalte
    std::cout << std::endl;
  }  // i - Zeile
  std::cout << std::endl;
}

int JacobiSIMD(__m128* uM128, const float vx, const float vy, const int nDiv4) {
  const int n = nDiv4 * 4;  // n immer durch 4 teilbar
  const float h = 1.0f / (n + 1.0f);
  const int maxNumIter = 2000;
  const float tol = 0.0000000001f;
  int Iter, i, j;  // Schleifenzähler
  float fourFloat[4];

  __m128* tempUM128;  // Vektor für u, zeilenweise Vektorisierung
  posix_memalign((void**)&tempUM128, 16, sizeof(float) * (n + 2) * (n + 8));
  __m128* uNewPtrM128 = tempUM128;  // um zwischen altem und neuem u zu wechseln
  __m128* uOldPtrM128 = uM128;
  float* uOld;  // um auf uM128/tempUM128 als float array zugreifen zu können

  // Alle u(x, y) = 0
  for (i = 0; i < n + 2; i++) {        // Zeilen
    for (j = 0; j < nDiv4 + 2; j++) {  // Spalten
      uM128[i * (nDiv4 + 2) + j] = _mm_setzero_ps();
      tempUM128[i * (nDiv4 + 2) + j] = _mm_set1_ps(1.0f);
    }  // j - Spalte
  }    // i - Zeile

  // Stern für Differenzen:
  const float a12 = -1.0f / h / h - vy / h;
  const float a21 = -1.0f / h / h - vx / h;
  const float a22 = 1.0f / (4.0f / h / h + vx / h + vy / h);  // Division incl.
  const float a23 = -1.0f / h / h;
  const float a32 = -1.0f / h / h;

  // Multiplikatoren für Stern in Vektorform:
  __m128 *a12M128, *a32M128, *a21M128, *a23M128, *a22M128;
  posix_memalign((void**)&a12M128, 16, sizeof(float) * 4);
  posix_memalign((void**)&a32M128, 16, sizeof(float) * 4);
  posix_memalign((void**)&a21M128, 16, sizeof(float) * 4);
  posix_memalign((void**)&a23M128, 16, sizeof(float) * 4);
  posix_memalign((void**)&a22M128, 16, sizeof(float) * 4);
  (*a12M128) = _mm_set1_ps(a12);
  (*a32M128) = _mm_set1_ps(a32);
  (*a21M128) = _mm_set1_ps(a21);
  (*a23M128) = _mm_set1_ps(a23);
  (*a22M128) = _mm_set1_ps(a22);

  // Vektoren für linke/rechte-Nachbarelemente (selbe Zeile in u)
  __m128 *leftUM128, *rightUM128;
  posix_memalign((void**)&leftUM128, 16, sizeof(float) * 4);
  posix_memalign((void**)&rightUM128, 16, sizeof(float) * 4);

  // Vektoren für Produkte u*Stern
  __m128 *leftProdM128, *rightProdM128, *upperProdM128, *lowerProdM128;
  posix_memalign((void**)&leftProdM128, 16, sizeof(float) * 4);
  posix_memalign((void**)&rightProdM128, 16, sizeof(float) * 4);
  posix_memalign((void**)&upperProdM128, 16, sizeof(float) * 4);
  posix_memalign((void**)&lowerProdM128, 16, sizeof(float) * 4);

  // Vektor für f=(1,1,1,1), f=(0,0,0,0), maxFehler, tol und temp:
  __m128 *oneM128, *zeroM128, *tolM128, *maxFehlerM128, *tempM128;
  posix_memalign((void**)&oneM128, 16, sizeof(float) * 4);
  (*oneM128) = _mm_set1_ps(1.0f);
  posix_memalign((void**)&zeroM128, 16, sizeof(float) * 4);
  (*zeroM128) = _mm_setzero_ps();
  posix_memalign((void**)&tolM128, 16, sizeof(float) * 4);
  (*tolM128) = _mm_set1_ps(tol);
  posix_memalign((void**)&maxFehlerM128, 16, sizeof(float) * 4);
  posix_memalign((void**)&tempM128, 16, sizeof(float) * 4);

  for (Iter = 0; Iter < maxNumIter; Iter++) {
    (*maxFehlerM128) = _mm_setzero_ps();  // Abbruchkriterium
    for (i = 1; i < n + 1; i++) {
      for (j = 1; j < nDiv4 + 1; j++) {  // zeilenweise vektorisiert

        uOld = (float*)uOldPtrM128;  // zum Heraussuchen linke/rechte Nachbarn
        (*leftUM128) = _mm_set_ps(
            uOld[i * (n + 8) + j * 4 + 2], uOld[i * (n + 8) + j * 4 + 1],
            uOld[i * (n + 8) + j * 4], uOld[i * (n + 8) + j * 4 - 1]);
        (*rightUM128) = _mm_set_ps(
            uOld[i * (n + 8) + j * 4 + 4], uOld[i * (n + 8) + j * 4 + 3],
            uOld[i * (n + 8) + j * 4 + 2], uOld[i * (n + 8) + j * 4 + 1]);

        (*leftProdM128) = _mm_mul_ps((*leftUM128), (*a21M128));
        (*rightProdM128) = _mm_mul_ps((*rightUM128), (*a23M128));
        (*upperProdM128) =
            _mm_mul_ps(uOldPtrM128[(i - 1) * (nDiv4 + 2) + j], (*a12M128));
        (*lowerProdM128) =
            _mm_mul_ps(uOldPtrM128[(i + 1) * (nDiv4 + 2) + j], (*a32M128));

        // Iterationsschritt:
        (*tempM128) = _mm_sub_ps((*oneM128), (*upperProdM128));
        (*tempM128) = _mm_sub_ps((*tempM128), (*lowerProdM128));
        (*tempM128) = _mm_sub_ps((*tempM128), (*leftProdM128));
        (*tempM128) = _mm_sub_ps((*tempM128), (*rightProdM128));
        uNewPtrM128[i * (nDiv4 + 2) + j] = _mm_mul_ps((*tempM128), (*a22M128));

        // Maximalwert (Absolutwert (Fehler)):
        (*tempM128) = _mm_sub_ps(uNewPtrM128[i * (nDiv4 + 2) + j],
                                 uOldPtrM128[i * (nDiv4 + 2) + j]);
        (*maxFehlerM128) = _mm_max_ps((*maxFehlerM128), (*tempM128));  // abs..
        (*tempM128) = _mm_sub_ps((*zeroM128), (*tempM128));            // abs..
        (*maxFehlerM128) = _mm_max_ps((*maxFehlerM128), (*tempM128));  // abs..
      }  // i - Zeilen von u
    }    // j - Spalten von u
    // Abbruchkriterium überprüfen:
    (*tempM128) = _mm_cmpgt_ps((*maxFehlerM128), (*tolM128));
    (*tempM128) = _mm_hadd_ps((*tempM128), (*tempM128));
    (*tempM128) = _mm_hadd_ps((*tempM128), (*tempM128));
    if (((float*)tempM128)[0] == 0.0f) return Iter;  // Genauigkeit erreicht
    std::swap(uOldPtrM128, uNewPtrM128);  // Vertauschen der Vektoren für u
  }                                       // Iter - Iterationen
  return Iter;                            // max. Anzahl an Iterationen erreicht
}

int main(int argc, char const* argv[]) {
  int nDiv4 = 256;  // Anzahl innere Gitterpunkte /4 (eine Richtung)
  int v0 = 50;

  // ggf. übergebene Parameter:
  if (argc > 2) {
    nDiv4 = atoi(argv[1]);
    v0 = atoi(argv[2]);
  }

  const int n = nDiv4 * 4;  // n immer durch 4 teilbar
  const float h = 1.0f / (n + 1.0f);
  const float vx = v0 / sqrt(5.0f);
  const float vy = 2.0f * v0 / sqrt(5.0f);
  __m128* uM128;  // vektor für u, zeilenweise vektorisiert
  posix_memalign((void**)&uM128, 16, sizeof(float) * (n + 2) * (n + 8));

  double startTime = clock();
  int ret = JacobiSIMD(uM128, vx, vy, nDiv4);
  double endTime = clock();

  // Ausgabe des Ergebnisses :
  // showMatrix(uM128, n + 2, n + 8);  // z.B. > u.dat
  std::cout << "Anzahl der Itterationen: " << ret << std::endl;
  std::cout << "Berechnungsdauer: " << endTime - startTime << std::endl;

  return 1;
}
