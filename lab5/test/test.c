#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fenv.h>
#include <float.h>

#define N 1024
int unrolling_factor = 8;
int tiling_factor = 64;

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

void multiply_double_sumator_kahana(double **d, double **e, long double **g) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            long double suma = 0.0;
            long double error = 0.0; // Błąd akumulacji
            for (int k = 0; k < N; k++) {
                long double y = d[i][k] * e[k][j] - error; // Obliczanie iloczynu z uwzględnieniem błędu
                long double t = suma + y;
                error = (t - suma) - y; // Obliczanie nowego błędu
                suma = t; // Aktualizacja sumy
            }
            g[i][j] = suma; // Zapisz wynik
        }
    }
}

void quicksort_long_double(long double *array, int low, int high) {
    if (low < high) {
        long double pivot = array[(low + high) / 2];
        int i = low, j = high;
        while (i <= j) {
            while (fabsl(array[i]) < fabsl(pivot)) i++;
            while (fabsl(array[j]) > fabsl(pivot)) j--;
            if (i <= j) {
                long double temp = array[i];
                array[i] = array[j];
                array[j] = temp;
                i++;
                j--;
            }
        }
        if (low < j) quicksort_long_double(array, low, j);
        if (i < high) quicksort_long_double(array, i, high);
    }
}

void multiply_double_sorted_sumator_kahana(double **d, double **e, long double **g) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            long double akumulator[N]; // Tablica do przechowywania iloczynów
            for (int k = 0; k < N; k++) {
                akumulator[k] = d[i][k] * e[k][j];
            }
            // Sortowanie iloczynów
            quicksort_long_double(akumulator, 0, N - 1);
            
            // Użycie algorytmu Kahana do sumowania
            long double suma = 0.0;
            long double error = 0.0; // Błąd akumulacji
            for (int k = 0; k < N; k++) {
                long double y = akumulator[k] - error; // Obliczanie iloczynu z uwzględnieniem błędu
                long double t = suma + y;
                error = (t - suma) - y; // Obliczanie nowego błędu
                suma = t; // Aktualizacja sumy
            }
            g[i][j] = suma; // Zapisz wynik
        }
    }
}

void calculate_error_long(long double **reference, long double **result, 
                     long double *min_error_long, long double *max_error_long, long double *mean_error_long) {
    *max_error_long = 0.0;
    *min_error_long = LDBL_MAX;// Ustawienie na maksymalną wartość dla double
    long double total_error = 0.0;
    int count = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Sprawdzenie, czy wartość referencyjna jest różna od zera, aby uniknąć dzielenia przez zero
            if (fabsl(reference[i][j]) > 0) {
                long double relative_error = fabsl(result[i][j] - reference[i][j]) / fabsl(reference[i][j]);
                total_error += relative_error;

                // Aktualizacja maksymalnego błędu względnego
                if (relative_error > *max_error_long) {
                    *max_error_long = relative_error;
                }

                // Aktualizacja minimalnego błędu względnego
                if (relative_error < *min_error_long) {
                    *min_error_long = relative_error;
                }

                count++;
            }
        }
    }

    // Obliczanie średniego błędu względnego
    if (count > 0) {
        *mean_error_long = total_error / count;
    } else {
        *mean_error_long = 0.0; // Brak wartości do obliczenia średniego błędu
        *min_error_long = 0.0; // Ustawienie na 0, gdy brak wartości
        *max_error_long = 0.0; // Ustawienie na 0, gdy brak wartości
    }
}

void reset(int **c, long double **g) {      // Reset macierze wynikowe
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c[i][j] = 0;
            g[i][j] = 0;
        }
    }
}

int main() {
    double range = 1000;
    double MIN_d = -range, MAX_d = range;
    unsigned long long start, end;
    unsigned int lo, hi;
    srand(time(NULL));

    long double min_error_long, max_error_long, mean_error_long;
    
    unsigned long long time_int = 0, time_double = 0;
    long double min_error_sum = 0, max_error_sum = 0, mean_error_sum = 0;
    
    int **a = (int **)malloc(N * sizeof(int *));
    int **b = (int **)malloc(N * sizeof(int *));
    int **c = (int **)malloc(N * sizeof(int *));
    int **b_transposed = malloc(N * sizeof(int *));
    double **d = malloc(N * sizeof(double*));
    double **e = malloc(N * sizeof(double*));
    
    long double **g = malloc(N * sizeof(long double*));
    long double **reference = malloc(N * sizeof(long double*));
    
    for (int i = 0; i < N; i++) {
        a[i] = (int *)malloc(N * sizeof(int));
        b[i] = (int *)malloc(N * sizeof(int));
        c[i] = (int *)malloc(N * sizeof(int));
        b_transposed[i] = (int *)malloc(N * sizeof(int));
        
        d[i] = malloc(N * sizeof(double));
        e[i] = malloc(N * sizeof(double));

        g[i] = malloc(N * sizeof(long double));
        reference[i] = malloc(N * sizeof(long double));
    }
    
    for (int r = 0; r < 10; r++) {
    // Inicjalizacja macierzy a i b oraz kopiowanie ich do d i e
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            a[i][j] = rand() % 10;
            b[i][j] = rand() % 10;
            c[i][j] = 0; // Inicjalizacja macierzy c

            d[i][j] = MIN_d + (MAX_d - MIN_d) * ((double)rand())/RAND_MAX;
            e[i][j] = MIN_d + (MAX_d - MIN_d) * ((double)rand())/RAND_MAX;
            
            g[i][j] = 0; // Inicjalizacja macierzy g
        }
    }
    
    multiply_double_sorted_sumator_kahana(d, e, reference); // Obliczanie macierzy referencyjnej do obliczania błędu
    printf("Macierz referencyjna obliczona\n\n");
    
    // Pomiar czasu dla mnożenia z optymalizacją cache, rozwinięciem pętli i blokowaniem
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_cache_unrolled_blocked(a, b, c, b_transposed, tiling_factor);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia z wszystkimi optymalizacjami (B=%d): %llu cykli\n\n", tiling_factor, end - start);
    time_int += (end - start);
 
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_double_sumator_kahana(d, e, g);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia double_sorted_sumator_kahana: %llu cykli\n", end - start);
    
    calculate_error_long(reference, g, &min_error_long, &max_error_long, &mean_error_long);
    printf("Minimalny błąd względny: %Le\n", min_error_long);
    printf("Maksymalny błąd względny: %Le\n", max_error_long);
    printf("Średni błąd względny: %Le\n\n", mean_error_long);
    
    time_double += (end - start);
    min_error_sum += min_error_long;
    max_error_sum += max_error_long;
    mean_error_sum += mean_error_long;
    reset(c, g);
    }
    
    printf("Średni czas mnożenia z wszystkimi optymalizacjami (B=%d): %llu cykli\n\n", tiling_factor, (time_int/100));
    printf("Średni czas mnożenia double_sorted_sumator_kahana: %llu cykli\n", (time_double/100));
    printf("Średni minimalny błąd względny: %Le\n", (min_error_sum/100));
    printf("Średni maksymalny błąd względny: %Le\n", (max_error_sum/100));
    printf("Średni średni błąd względny: %Le\n\n", (mean_error_sum/100));

    for (int i = 0; i < N; i++) {
        free(a[i]); free(b[i]); free(c[i]);
        free(b_transposed[i]);
        free(d[i]); free(e[i]);
        free(g[i]); free(reference[i]);
    }
    free(a); free(b); free(c);
    free(b_transposed);
    free(d); free(e);
    free(g); free(reference);

    return 0;
} 
