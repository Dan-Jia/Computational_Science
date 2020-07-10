#include <algorithm>
#include <cmath>
#include <iostream>

void showMatrix(float* M, int n, int m) {  // zum Auslese nach Matlab z.B.
  std::cout << std::fixed;
  std::cout << std::endl;
  for (int i = 0; i < n; i++) {    // Zeilen
    for (int j = 0; j < m; j++) {  // Spalten
      std::cout << M[i * m + j] << " ";
    }  // j - Spalte
    std::cout << std::endl;
  }  // i - Zeile
  std::cout << std::endl;
}

int JacobiSkalar(float* u, float* tempU, const float vx, const float vy,
                 const int nDiv4) {
  const long int n = nDiv4 * 4;  // innere Gitterpunkte
  const float h = 1.0f / (n + 1.0f);
  const int maxNumIter = 2000;
  float* uNewPtr = tempU;
  float* uOldPtr = u;
  float maxFehler;      // für Abbruchkriterium
  long int Iter, i, j;  // Schleifenzähler

  // Alle u(x,y)=0;
  for (i = 0; i < n + 2; i++)
    for (j = 0; j < n + 2; j++) {
      tempU[i * (n + 2) + j] = 1.0f;
      u[i * (n + 2) + j] = 0.0f;
    }

  // Stern für linksseitige Differenzen:
  const float a12 = -1.0f / h / h - vy / h;
  const float a21 = -1.0f / h / h - vx / h;
  const float a22 = 1 / (4.0f / h / h + vx / h + vy / h);  // incl. 1/..
  const float a23 = -1.0f / h / h;
  const float a32 = -1.0f / h / h;

  // Jacobi Iteration:
  for (Iter = 0; Iter < maxNumIter; Iter++) {
    maxFehler = 0.0f;
    for (i = 1; i < n + 1; i++) {    // Zeilen
      for (j = 1; j < n + 1; j++) {  // Spalten
        uNewPtr[i * (n + 2) + j] =
            (1.0f - uOldPtr[i * (n + 2) + (j - 1)] * a21 -
             uOldPtr[i * (n + 2) + (j + 1)] * a23 -
             uOldPtr[(i - 1) * (n + 2) + j] * a12 -
             uOldPtr[(i + 1) * (n + 2) + j] * a32) *
            a22;
        maxFehler =
            fmaxf(abs(uNewPtr[i * n + j] - uOldPtr[i * n + j]), maxFehler);
      }                                     // i - Zeilen von u
    }                                       // j - Spalten von u
    if (maxFehler < 0.00001f) return Iter;  // Genauigkeit erreicht
    std::swap(uOldPtr, uNewPtr);
  }  // Iter

  return -999;  // max. Anzahl an Iterationen erreicht
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
  float* u = new float[(n + 2) * (n + 2)];
  float* uTemp = new float[(n + 2) * (n + 2)];
  const float h = 1.0f / (n + 1.0f);
  const float vx = v0 / sqrt(5.0f);
  const float vy = 2.0f * v0 / sqrt(5.0f);

  double startTime = clock();
  int ret = JacobiSkalar(u, uTemp, vx, vy, nDiv4);
  double endTime = clock();

  // Ausgabe des Ergebnisses :
  // showMatrix(u, n + 2, n + 2);  // z.B. > u.dat
  std::cout << "Anzahl der Itterationen: " << ret << std::endl;
  std::cout << "Berechnungsdauer: " << endTime - startTime << std::endl;

  return 1;
}
