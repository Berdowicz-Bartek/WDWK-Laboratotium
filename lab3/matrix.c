#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 64
int a[N][N], b[N][N], c[N][N];

void multiply_naive(int a[N][N], int b[N][N], int c[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c[i][j] = 0; // Inicjalizacja elementu
            for (int k = 0; k < N; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void multiply_cache(int a[N][N], int b[N][N], int c[N][N]) {
    // Transponowanie macierzy b
    int b_transposed[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b_transposed[j][i] = b[i][j];
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c[i][j] = 0; // Inicjalizacja elementu
            for (int k = 0; k < N; k++) {
                c[i][j] += a[i][k] * b_transposed[j][k];
            }
        }
    }
}

void multiply_unrolled(int a[N][N], int b[N][N], int c[N][N]) {
    int unrolling_factor = 2;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c[i][j] = 0; // Inicjalizacja elementu
            for (int k = 0; k < N; k += unrolling_factor) {
                c[i][j] += a[i][k] * b[k][j];
                c[i][j] += a[i][k + 1] * b[k + 1][j];
            }
        }
    }
}

void multiply_blocked(int a[N][N], int b[N][N], int c[N][N], int B) {
    for (int i = 0; i < N; i += B) {
        for (int j = 0; j < N; j += B) {
            for (int k = 0; k < N; k += B) {
                for (int ii = i; ii < i + B; ii++) {
                    for (int jj = j; jj < j + B; jj++) {
                        for (int kk = k; kk < k + B; kk++) {
                            c[ii][jj] += a[ii][kk] * b[kk][jj];
                        }
                    }
                }
            }
        }
    }
}

void multiply_cache_unrolled(int a[N][N], int b[N][N], int c[N][N]) {
    // Transponowanie macierzy b
    int b_transposed[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b_transposed[j][i] = b[i][j];
        }
    }

    int unrolling_factor = 2;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c[i][j] = 0; // Inicjalizacja elementu
            for (int k = 0; k < N; k += unrolling_factor) {
                c[i][j] += a[i][k] * b_transposed[j][k];
                c[i][j] += a[i][k + 1] * b_transposed[j][k + 1];
            }
        }
    }
}

void multiply_cache_blocked(int a[N][N], int b[N][N], int c[N][N], int B) {
    // Transponowanie macierzy b
    int b_transposed[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b_transposed[j][i] = b[i][j];
        }
    }

    for (int i = 0; i < N; i += B) {
        for (int j = 0; j < N; j += B) {
            for (int k = 0; k < N; k += B) {
                for (int ii = i; ii < i + B; ii++) {
                    for (int jj = j; jj < j + B; jj++) {
                        for (int kk = k; kk < k + B; kk++) {
                            c[ii][jj] += a[ii][kk] * b_transposed[jj][kk];
                        }
                    }
                }
            }
        }
    }
}

void multiply_unrolled_blocked(int a[N][N], int b[N][N], int c[N][N], int B) {
    for (int i = 0; i < N; i += B) {
        for (int j = 0; j < N; j += B) {
            for (int k = 0; k < N; k += B) {
                for (int ii = i; ii < i + B; ii++) {
                    for (int jj = j; jj < j + B; jj++) {
                        for (int kk = k; kk < k + B; kk += 2) {
                            c[ii][jj] += a[ii][kk] * b[kk][jj];
                            c[ii][jj] += a[ii][kk + 1] * b[kk + 1][jj];
                        }
                    }
                }
            }
        }
    }
}

void multiply_cache_unrolled_blocked(int a[N][N], int b[N][N], int c[N][N], int B) {
    // Transponowanie macierzy b
    int b_transposed[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b_transposed[j][i] = b[i][j];
        }
    }

    for (int i = 0; i < N; i += B) {
        for (int j = 0; j < N; j += B) {
            for (int k = 0; k < N; k += B) {
                for (int ii = i; ii < i + B; ii++) {
                    for (int jj = j; jj < j + B; jj++) {
                        for (int kk = k; kk < k + B; kk += 2) {
                            c[ii][jj] += a[ii][kk] * b_transposed[jj][kk];
                            c[ii][jj] += a[ii][kk + 1] * b_transposed[jj][kk + 1];
                        }
                    }
                }
            }
        }
    }
}

void print_matrix(int matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", matrix[i][j]);
        }
        printf("\n");
    }
}

void reset() {      // Reset macierzy c
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c[i][j] = 0;
        }
    }
}

int main() {
    unsigned long long start, end, time;
    unsigned int lo, hi;
    int B = 4;   // Przykładowy tiling factor
    
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
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_naive(a, b, c);
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia naiwnym: %llu cykli\n", end - start);
    
    reset();

    // Pomiar czasu dla mnożenia z optymalizacją cache
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_cache(a, b, c);
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z optymalizacją cache: %llu cykli\n", end - start);

    reset();

    // Pomiar czasu dla mnożenia z rozwinięciem pętli
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_unrolled(a, b, c);
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z rozwinięciem pętli: %llu cykli\n", end - start);

    reset();

    // Pomiar czasu dla mnożenia z blokowaniem
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_blocked(a, b, c, B);
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z blokowaniem (B=%d): %llu cykli\n", B, end - start);
    
    reset();
    printf("\n");

    // Pomiar czasu dla mnożenia z optymalizacją cache i rozwinięciem pętli
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_cache_unrolled(a, b, c);
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z optymalizacją cache i rozwinięciem pętli: %llu cykli\n", end - start);
    
    reset();
    
    // Pomiar czasu dla mnożenia z optymalizacją cache i blokowaniem
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_cache_blocked(a, b, c, B);
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z optymalizacją cache i blokowaniem (B=%d): %llu cykli\n", B, end - start);
    
    reset();
    
    // Pomiar czasu dla mnożenia z rozwinięciem pętli i blokowaniem
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_unrolled_blocked(a, b, c, B);
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z rozwinięciem pętli i blokowaniem (B=%d): %llu cykli\n", B, end - start);
    
    reset();
    
    // Pomiar czasu dla mnożenia z optymalizacją cache, rozwinięciem pętli i blokowaniem
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_cache_unrolled_blocked(a, b, c, B);
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z wszystkimi optymalizacjami (B=%d): %llu cykli\n", B, end - start);
    
    return 0;
} 
