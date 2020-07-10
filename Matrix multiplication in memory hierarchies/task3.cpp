#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

using base_t = std::int64_t;

void run(const base_t* buffer, int size, int stride) {
  volatile base_t sum = 0;
  for (int i = 0; i < size; ++i) sum += buffer[i * stride];
}

int main(int argc, char const* argv[]) {
  int MIN_STRIDE = 1;
  int MAX_STRIDE = 16;
  int STRIDE_STEP = 1;

  int MIN_SIZE = 2048;
  int MAX_SIZE = 67108864;
  int SIZE_FACT = 2;

  int NUM_RUNS_FACT = MAX_SIZE;

  std::vector<std::tuple<int, int, double>> measurements;

  for (int stride = MIN_STRIDE; stride <= MAX_STRIDE; stride += STRIDE_STEP) {
    for (int size = MIN_SIZE; size <= MAX_SIZE; size *= SIZE_FACT) {
      const int num_elements = size / sizeof(base_t);
      std::vector<base_t> buffer(num_elements * stride, 42);
      using namespace std::chrono;
      auto start = high_resolution_clock::now();
      const int NUM_RUNS = NUM_RUNS_FACT / size;
      for (int i = 0; i < NUM_RUNS; ++i)
        run(buffer.data(), num_elements, stride);
      auto end = high_resolution_clock::now();
      const double elapsed_seconds = duration<double>(end - start).count();
      const double bandwidth = 1e-6 * size / (elapsed_seconds / NUM_RUNS);
      measurements.emplace_back(stride, size, bandwidth);
    }
  }

  std::ofstream measurements_file("measurements.dat");
  for (const auto& e : measurements)
    measurements_file << std::get<0>(e) << ' ' << std::get<1>(e) << ' '
                      << std::get<2>(e) << std::endl;
}
