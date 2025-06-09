#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define KB(x) ((x) * 1024)
#define MB(x) ((x) * 1024 * 1024)

#define N_REPEATS 5
#define CPU_FREQ_HZ 2592000000.0  // 2.592 GHz

size_t sizes[] = {
    KB(1), KB(4), KB(16), KB(64), KB(256),
    MB(1), MB(4), MB(16), MB(64), MB(128), MB(256)
};

size_t strides[] = { 4, 16, 64, 256, 512, 1024, 2048 };

volatile uint64_t sum;

// Funkcja do pomiaru czasu przy pomocy rdtsc
static inline uint64_t rdtsc_start() {
    unsigned int lo, hi;
    __asm__ volatile (
        "cpuid\n\t"
        "rdtsc\n\t"
        : "=a"(lo), "=d"(hi)
        :
        : "%rbx", "%rcx"
    );
    return ((uint64_t)hi << 32) | lo;
}

static inline uint64_t rdtsc_end() {
    unsigned int lo, hi;
    __asm__ volatile (
        "rdtscp\n\t"
        : "=a"(lo), "=d"(hi)
        :
        : "%rcx"
    );
    __asm__ volatile ("cpuid\n\t" : : : "%rax", "%rbx", "%rcx", "%rdx");
    return ((uint64_t)hi << 32) | lo;
}

void run_test(size_t size, size_t stride, int cold_cache) {
    char *array = (char *)malloc(size);
    if (!array) {
        fprintf(stderr, "Allocation failed for size %zu\n", size);
        return;
    }

    for (size_t i = 0; i < size; i++) {
        array[i] = (char)(i % 256);
    }

    uint64_t best_cycles = UINT64_MAX;
    uint64_t temp = 0;

    for (int repeat = 0; repeat < N_REPEATS; repeat++) {
        if (cold_cache) {
            memset(array, repeat, size);
        }

        uint64_t start = rdtsc_start();

        uint64_t local_sum = 0;
        for (size_t i = 0; i < size; i += stride) {
            local_sum += array[i];
        }

        uint64_t end = rdtsc_end();
        if (end - start < best_cycles) {
            best_cycles = end - start;
            temp = local_sum;
        }
    }

    sum = temp; // zapobiega optymalizacji

    // Obliczanie przepustowości
    size_t accesses = size / stride;
    size_t bytes_accessed = accesses * sizeof(char);
    double time_seconds = best_cycles / CPU_FREQ_HZ;
    double bandwidth_gbps = (bytes_accessed / time_seconds) / (1024.0 * 1024.0 * 1024.0);

    printf("%zu,%zu,%s,%lu,%.6f\n",
           size, stride, cold_cache ? "cold" : "warm", best_cycles, bandwidth_gbps);

    free(array);
}

int main() {
    printf("Size(Bytes),Stride(Bytes),CacheState,Cycles,Bandwidth_GBps\n");

    for (size_t s = 0; s < sizeof(sizes)/sizeof(sizes[0]); s++) {
        for (size_t st = 0; st < sizeof(strides)/sizeof(strides[0]); st++) {
            run_test(sizes[s], strides[st], 1); // Cold cache
            run_test(sizes[s], strides[st], 0); // Warm cache
        }
    }

    return 0;
}

