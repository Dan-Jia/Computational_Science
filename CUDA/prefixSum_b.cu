#include <cuda.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
using std::setw;

const int N = 5;
const int threadsPerBlock = N;
int blocksPerGrid = 1;

// device code
__global__ void prefixSum(float* x, float* c) {
  __shared__ float
      cache[2 * threadsPerBlock];  // declaring a array in shared memory
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int cacheIdx = threadIdx.x;  // cacheIdx equal the threadIdx in every block

  if (tid < N) {
    cache[cacheIdx] = x[tid];
  }

  // reduction step: the code below performs iterative scan on cache
  for (int stride = 1; stride <= threadsPerBlock; stride *= 2) {
    __syncthreads();
    int index = (threadIdx.x + 1) * stride * 2 - 1;
    if (index < 2 * threadsPerBlock) {
      cache[index] +=
          cache[index - stride];  // index is alway bigger than stride
    }
    __syncthreads();
  }
  // threadIdx.x+1 = 1,2,3,4....
  // stride index = 1,3,5,7...

  for (int stride = threadsPerBlock / 2; stride > 0; stride /= 2) {
    __syncthreads();
    int index = (threadIdx.x + 1) * stride * 2 - 1;
    if (index < 2 * threadsPerBlock) {
      cache[index] +=
          cache[index - stride];  // index is alway bigger than stride
    }
    __syncthreads();
  }

  // reduction reverse phase
  int stride = 1;
  while (stride <= threadsPerBlock) {
    int index = (threadIdx.x + 1) * stride * 2 - 1;
    if (index < 2 * threadsPerBlock) {
      cache[index + stride] += cache[index];
    }
  }
  __syncthreads();

  if (tid < N) {
    c[tid] = cache[threadIdx.x];
  }
}

// host code
int main() {
  size_t size = N * sizeof(float);
  size_t size_c = blocksPerGrid * sizeof(float);

  // Allocate input vectors h_a and h_b in host memory
  float* h_x = (float*)malloc(size);
  float* h_c = (float*)malloc(size_c);

  // Initialize input vectors
  for (int i = 0; i < N; ++i) {
    h_x[i] = i;
  }

  // Allocate vectors in device memory
  float* d_x;
  float* d_c;

  cudaMalloc(&d_x, size);
  cudaMalloc(&d_c, size_c);

  // copy vectors from host memory to device memory
  cudaMemcpy(d_x, h_x, size, cudaMemcpyHostToDevice);

  // Invoke kernel
  prefixSum<<<blocksPerGrid, threadsPerBlock>>>(d_x, d_c);

  // copy result from device memory to host memory
  // h_c contains the result in host memory
  cudaMemcpy(h_c, d_c, size_c, cudaMemcpyDeviceToHost);

  // output the result
  std::cout << "Element" << setw(13) << "Value" << std::endl;
  for (int i = 0; i < N; ++i) {
    std::cout << setw(7) << i << setw(13) << h_c[i] << std::endl;
  }

  // Free device memory
  cudaFree(d_x);
  cudaFree(d_c);

  // Free host memory
  free(h_x);
  free(h_c);

  cudaThreadExit();
  return 0;
}
