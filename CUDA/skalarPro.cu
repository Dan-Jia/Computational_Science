#include <cuda.h>
#include <stdio.h>

const int N = 524288;
const int threadsPerBlock = 1024;
int blocksPerGrid = (N + threadsPerBlock) / threadsPerBlock - 1;

// device code
__global__ void skalarPro(float* a, float* b, float* c) {
  __shared__ float cache[threadsPerBlock];  // declaring a array for every block
                                            // in shared memory

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int cacheIdx = threadIdx.x;  // cacheIdx equal the threadIdx in every block

  // save the result of every thread in cache[]
  // if the computing can not once finish
  float temp = 0;
  while (tid < N) {
    temp += a[tid] * b[tid];
    tid += blockDim.x * gridDim.x;  // adding temp per all threads
  }

  cache[cacheIdx] = temp;  // save the result of every thread in cache[cacheIdx]
  __syncthreads();

  // reduce the sum in every cache to cache[0]
  int i = blockDim.x / 2;
  while (i != 0) {
    if (cacheIdx < i) {
      cache[cacheIdx] += cache[cacheIdx + i];
    }
    __syncthreads();
    i /= 2;
  }

  // store the result of cache[0] into global variable c[]
  if (cacheIdx == 0) {
    c[blockIdx.x] = cache[0];
  }
}

// host code
int main() {
  size_t size = N * sizeof(float);
  size_t size_c = blocksPerGrid * sizeof(float);

  // Allocate input vectors h_a and h_b in host memory
  float* h_a = (float*)malloc(size);
  float* h_b = (float*)malloc(size);
  float* h_c = (float*)malloc(size_c);

  // Initialize input vectors
  for (int i = 0; i < N; i++) {
    h_a[i] = i + 1;
    h_b[i] = N - i;
  }

  // Allocate vectors in device memory
  float* d_a;
  float* d_b;
  float* d_c;

  cudaMalloc(&d_a, size);
  cudaMalloc(&d_b, size);
  cudaMalloc(&d_c, size_c);

  // copy vectors from host memory to device memory
  cudaMemcpy(d_a, h_a, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, h_b, size, cudaMemcpyHostToDevice);

  // Invoke kernel
  skalarPro<<<blocksPerGrid, threadsPerBlock>>>(d_a, d_b, d_c);

  // copy result from device memory to host memory
  // h_c contains the result in host memory
  cudaMemcpy(h_c, d_c, size_c, cudaMemcpyDeviceToHost);

  // reduce the sum in every h_c[i] to cache[0]
  int i = blocksPerGrid / 2;
  while (i != 0) {
    for (int j = 0; j < blocksPerGrid; j++) {
      if (j < i) {
        h_c[j] += h_c[j + i];
      }
    }
    i /= 2;
  }

  printf("Das Ergebnis vom Skalarprodukt: %f\n", h_c[0]);

  // Free device memory
  cudaFree(d_a);
  cudaFree(d_b);
  cudaFree(d_c);

  // Free host memory
  free(h_a);
  free(h_b);
  free(h_c);

  cudaThreadExit();
  return 0;
}
