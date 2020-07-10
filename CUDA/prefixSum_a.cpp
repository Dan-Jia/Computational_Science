#include <stdlib.h>
#include <iomanip>
#include <iostream>
using std::setw;

void prefixSum(float* x, float* s, int N) {
  s[0] = x[0];
  for (int i = 1; i < N; ++i) {
    s[i] = s[i - 1] + x[i];
  }

  // output all values of the array
  std::cout << "Element" << setw(13) << "Value" << std::endl;
  for (int i = 0; i < N; ++i) {
    std::cout << setw(7) << i << setw(13) << s[i] << std::endl;
  }
}

int main() {
  int N = 5;
  float* x = (float*)malloc(N * sizeof(float));
  float* s = (float*)malloc(N * sizeof(float));
  for (int i = 0; i < N; ++i) {
    x[i] = i;
  }

  // Invoke function prefixSum
  prefixSum(x, s, N);

  free(x);
  free(s);

  return 0;
}
