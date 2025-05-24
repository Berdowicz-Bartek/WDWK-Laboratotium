#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 1024  //rozmiar macierzy
int unrolling_factor = 8;
int tiling_factor = 64;


void multiply_naive(int **a, int **b, int **c) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void multiply_cache(int **a, int **b, int **c, int **b_transposed) {
    // Transponowanie macierzy b
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b_transposed[j][i] = b[i][j];
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                c[i][j] += a[i][k] * b_transposed[j][k];
            }
        }
    }
}

void multiply_unrolled(int **a, int **b, int **c) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int sum = 0;
            for (int k = 0; k < N; k += unrolling_factor) {
                sum += a[i][k] * b[k][j];
                sum += a[i][k + 1] * b[k + 1][j];
                sum += a[i][k + 2] * b[k + 2][j];
                sum += a[i][k + 3] * b[k + 3][j];
                sum += a[i][k + 4] * b[k + 4][j];
                sum += a[i][k + 5] * b[k + 5][j];
                sum += a[i][k + 6] * b[k + 6][j];
                sum += a[i][k + 7] * b[k + 7][j];
            }
            c[i][j] = sum;
        }
    }
}

void multiply_blocked(int **a, int **b, int **c, int tiling_factor) {
    for (int i = 0; i < N; i += tiling_factor) {
        for (int j = 0; j < N; j += tiling_factor) {
            for (int k = 0; k < N; k += tiling_factor) {
                for (int ii = i; ii < i + tiling_factor; ii++) {
                    for (int jj = j; jj < j + tiling_factor; jj++) {
                        for (int kk = k; kk < k + tiling_factor; kk++) {
                            c[ii][jj] += a[ii][kk] * b[kk][jj];
                        }
                    }
                }
            }
        }
    }
}

void multiply_cache_unrolled(int **a, int **b, int **c, int **b_transposed) {
    // Transponowanie macierzy b
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b_transposed[j][i] = b[i][j];
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int sum = 0;
            for (int k = 0; k < N; k += unrolling_factor) {
                sum += a[i][k] * b_transposed[j][k];
                sum += a[i][k + 1] * b_transposed[j][k + 1];
                sum += a[i][k + 2] * b_transposed[j][k + 2];
                sum += a[i][k + 3] * b_transposed[j][k + 3];
                sum += a[i][k + 4] * b_transposed[j][k + 4];
                sum += a[i][k + 5] * b_transposed[j][k + 5];
                sum += a[i][k + 6] * b_transposed[j][k + 6];
                sum += a[i][k + 7] * b_transposed[j][k + 7];
            }
            c[i][j] = sum;
        }
    }
}

void multiply_cache_blocked(int **a, int **b, int **c, int **b_transposed, int tiling_factor) {
    // Transponowanie macierzy b
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b_transposed[j][i] = b[i][j];
        }
    }

    for (int i = 0; i < N; i += tiling_factor) {
        for (int j = 0; j < N; j += tiling_factor) {
            for (int k = 0; k < N; k += tiling_factor) {
                for (int ii = i; ii < i + tiling_factor; ii++) {
                    for (int jj = j; jj < j + tiling_factor; jj++) {
                        for (int kk = k; kk < k + tiling_factor; kk++) {
                            c[ii][jj] += a[ii][kk] * b_transposed[jj][kk];
                        }
                    }
                }
            }
        }
    }
}

void multiply_unrolled_blocked(int **a, int **b, int **c, int tiling_factor) {
    for (int i = 0; i < N; i += tiling_factor) {
        for (int j = 0; j < N; j += tiling_factor) {
            for (int k = 0; k < N; k += tiling_factor) {
                for (int ii = i; ii < i + tiling_factor; ii++) {
                    for (int jj = j; jj < j + tiling_factor; jj++) {
                        int sum = 0;
                        for (int kk = k; kk < k + tiling_factor; kk += unrolling_factor) {
                            sum += a[ii][kk] * b[kk][jj];
                            sum += a[ii][kk + 1] * b[kk + 1][jj];
                            sum += a[ii][kk + 2] * b[kk + 2][jj];
                            sum += a[ii][kk + 3] * b[kk + 3][jj];
                            sum += a[ii][kk + 4] * b[kk + 4][jj];
                            sum += a[ii][kk + 5] * b[kk + 5][jj];
                            sum += a[ii][kk + 6] * b[kk + 6][jj];
                            sum += a[ii][kk + 7] * b[kk + 7][jj];
                        }
                        c[ii][jj] = sum;
                    }
                }
            }
        }
    }
}

void multiply_cache_unrolled_blocked(int **a, int **b, int **c, int **b_transposed, int tiling_factor) {
    // Transponowanie macierzy b
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b_transposed[j][i] = b[i][j];
        }
    }

    for (int i = 0; i < N; i += tiling_factor) {
        for (int j = 0; j < N; j += tiling_factor) {
            for (int k = 0; k < N; k += tiling_factor) {
                for (int ii = i; ii < i + tiling_factor; ii++) {
                    for (int jj = j; jj < j + tiling_factor; jj++) {
                        int sum = 0;
                        for (int kk = k; kk < k + tiling_factor; kk += unrolling_factor) {
                            sum += a[ii][kk] * b_transposed[jj][kk];
                            sum += a[ii][kk + 1] * b_transposed[jj][kk + 1];
                            sum += a[ii][kk + 2] * b_transposed[jj][kk + 2];
                            sum += a[ii][kk + 3] * b_transposed[jj][kk + 3];
                            sum += a[ii][kk + 4] * b_transposed[jj][kk + 4];
                            sum += a[ii][kk + 5] * b_transposed[jj][kk + 5];
                            sum += a[ii][kk + 6] * b_transposed[jj][kk + 6];
                            sum += a[ii][kk + 7] * b_transposed[jj][kk + 7];
                        }
                        c[ii][jj] = sum;
                    }
                }
            }
        }
    }
}

void print_matrix(int **matrix) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

void reset(int **c) {      // Reset macierzy c
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c[i][j] = 0;
        }
    }
}

int main() {
    unsigned long long start, end;
    unsigned int lo, hi;
    
    // Dynamiczna alokacja pamięci dla macierzy a, b i c
    int **a = (int **)malloc(N * sizeof(int *));
    int **b = (int **)malloc(N * sizeof(int *));
    int **c = (int **)malloc(N * sizeof(int *));
    int **b_transposed = malloc(N * sizeof(int *));
    for (int i = 0; i < N; i++) {
        a[i] = (int *)malloc(N * sizeof(int));
        b[i] = (int *)malloc(N * sizeof(int));
        c[i] = (int *)malloc(N * sizeof(int));
        b_transposed[i] = (int *)malloc(N * sizeof(int));
    }
    
    // Inicjalizacja macierzy a i b
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            a[i][j] = rand() % 10;
            b[i][j] = rand() % 10;
            c[i][j] = 0; // Inicjalizacja macierzy c
        }
    }
    /*
    print_matrix(a);
    printf("\n");
    print_matrix(b);
    printf("\n");
    print_matrix(c);
*/
    // Pomiar czasu dla mnożenia naiwnym
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_naive(a, b, c);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia naiwnym: %llu cykli\n", end - start);

    reset(c);

    // Pomiar czasu dla mnożenia z optymalizacją cache
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_cache(a, b, c, b_transposed);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z optymalizacją cache: %llu cykli\n", end - start);

    reset(c);

    // Pomiar czasu dla mnożenia z rozwinięciem pętli
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_unrolled(a, b, c);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z rozwinięciem pętli: %llu cykli\n", end - start);

    reset(c);

    // Pomiar czasu dla mnożenia z blokowaniem
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_blocked(a, b, c, tiling_factor);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z blokowaniem (B=%d): %llu cykli\n", tiling_factor, end - start);
    
    reset(c);
    printf("\n");

    // Pomiar czasu dla mnożenia z optymalizacją cache i rozwinięciem pętli
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_cache_unrolled(a, b, c, b_transposed);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z optymalizacją cache i rozwinięciem pętli: %llu cykli\n", end - start);
    
    reset(c);
    
    // Pomiar czasu dla mnożenia z optymalizacją cache i blokowaniem
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_cache_blocked(a, b, c, b_transposed, tiling_factor);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z optymalizacją cache i blokowaniem (B=%d): %llu cykli\n", tiling_factor, end - start);
    
    reset(c);
    
    // Pomiar czasu dla mnożenia z rozwinięciem pętli i blokowaniem
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_unrolled_blocked(a, b, c, tiling_factor);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z rozwinięciem pętli i blokowaniem (B=%d): %llu cykli\n", tiling_factor, end - start);
    
    reset(c);
    
    // Pomiar czasu dla mnożenia z optymalizacją cache, rozwinięciem pętli i blokowaniem
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_cache_unrolled_blocked(a, b, c, b_transposed, tiling_factor);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z wszystkimi optymalizacjami (B=%d): %llu cykli\n", tiling_factor, end - start);
    
    // Zwolnienie pamięci
    for (int i = 0; i < N; i++) {
        free(a[i]);
        free(b[i]);
        free(c[i]);
        free(b_transposed[i]);
    }
    free(a);
    free(b);
    free(c);
    free(b_transposed);
    
    return 0;
} 
