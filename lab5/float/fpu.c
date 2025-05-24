#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fenv.h>
#include <float.h>

#define N 700

void multiply_float(float **a, float **b, float **c) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void multiply_double(double **d, double **e, double **f) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                f[i][j] += d[i][k] * e[k][j];
            }
        }
    }
}

void quicksort_float(float *array, int low, int high) {
    if (low < high) {
        float pivot = array[(low + high) / 2]; // Wybór pivotu jako mediany
        int i = low, j = high;
        while (i <= j) {
            while (fabsf(array[i]) < fabsf(pivot)) i++;
            while (fabsf(array[j]) > fabsf(pivot)) j--;
            if (i <= j) {
                float temp = array[i];
                array[i] = array[j];
                array[j] = temp;
                i++;
                j--;
            }
        }
        if (low < j) quicksort_float(array, low, j);
        if (i < high) quicksort_float(array, i, high);
    }
}

void quicksort_double(double *array, int low, int high) {
    if (low < high) {
        double pivot = array[(low + high) / 2];
        int i = low, j = high;
        while (i <= j) {
            while (fabs(array[i]) < fabs(pivot)) i++;
            while (fabs(array[j]) > fabs(pivot)) j--;
            if (i <= j) {
                double temp = array[i];
                array[i] = array[j];
                array[j] = temp;
                i++;
                j--;
            }
        }
        if (low < j) quicksort_double(array, low, j);
        if (i < high) quicksort_double(array, i, high);
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

void multiply_float_sorted(float **d, float **e, float **f) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float akumulator[N]; // Tablica do przechowywania iloczynów
            for (int k = 0; k < N; k++) {
                akumulator[k] = d[i][k] * e[k][j];
            }
            quicksort_float(akumulator, 0, N - 1);
            for (int k = 0; k < N; k++) {
                f[i][j] += akumulator[k];
            }
        }
    }
}

void multiply_double_sorted(double **d, double **e, double **f) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double akumulator[N]; // Tablica do przechowywania iloczynów
            for (int k = 0; k < N; k++) {
                akumulator[k] = d[i][k] * e[k][j];
            }
            quicksort_double(akumulator, 0, N-1);
            for (int k = 0; k < N; k++) {
                f[i][j] += akumulator[k];
            }
        }
    }
}

void multiply_float_sumator(float **a, float **b, double **f) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            f[i][j] = 0; // Inicjalizacja elementu
            for (int k = 0; k < N; k++) {
                f[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void multiply_double_sumator(double **d, double **e, long double **g) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            g[i][j] = 0; // Inicjalizacja elementu
            for (int k = 0; k < N; k++) {
                g[i][j] += d[i][k] * e[k][j];
            }
        }
    }
}

void multiply_float_kahana(float **a, float **b, float **c) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float suma = 0.0;
            float error = 0.0;   // Błąd akumulacji
            for (int k = 0; k < N; k++) {
                float y = a[i][k] * b[k][j] - error; // Obliczanie iloczynu z uwzględnieniem błędu
                float t = suma + y;
                error = (t - suma) - y; // Obliczanie nowego błędu
                suma = t; // Aktualizacja sumy
            }
            c[i][j] = suma;
        }
    }
}

void multiply_double_kahana(double **d, double **e, double **f) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double suma = 0.0;
            double error = 0.0;   // Błąd akumulacji
            for (int k = 0; k < N; k++) {
                double y = d[i][k] * e[k][j] - error; // Obliczanie iloczynu z uwzględnieniem błędu
                double t = suma + y;
                error = (t - suma) - y; // Obliczanie nowego błędu
                suma = t; // Aktualizacja sumy
            }
            f[i][j] = suma;
        }
    }
}

void multiply_float_sorted_sumator(float **a, float **b, double **f) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double akumulator[N]; // Tablica do przechowywania iloczynów
            for (int k = 0; k < N; k++) {
                akumulator[k] = a[i][k] * b[k][j];
            }
            quicksort_double(akumulator, 0, N-1); // Sortowanie akumulatora
            for (int k = 0; k < N; k++) {
                f[i][j] += akumulator[k]; // Sumowanie posortowanych wartości
            }
        }
    }
}

void multiply_double_sorted_sumator(double **d, double **e, long double **g) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            long double akumulator[N]; // Tablica do przechowywania iloczynów
            for (int k = 0; k < N; k++) {
                akumulator[k] = d[i][k] * e[k][j];
            }
            quicksort_long_double(akumulator, 0, N-1); // Sortowanie akumulatora
            for (int k = 0; k < N; k++) {
                g[i][j] += akumulator[k]; // Sumowanie posortowanych wartości
            }
        }
    }
}

void multiply_float_sumator_kahana(float **a, float **b, double **f) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double suma = 0.0;
            double error = 0.0; // Błąd akumulacji
            for (int k = 0; k < N; k++) {
                double y = a[i][k] * b[k][j] - error; // Obliczanie iloczynu z uwzględnieniem błędu
                double t = suma + y;
                error = (t - suma) - y; // Obliczanie nowego błędu
                suma = t; // Aktualizacja sumy
            }
            f[i][j] = suma; // Zapisz wynik
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

void multiply_float_sorted_kahana(float **a, float **b, float **c) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float akumulator[N]; // Tablica do przechowywania iloczynów
            // Mnożenie macierzy
            for (int k = 0; k < N; k++) {
                akumulator[k] = a[i][k] * b[k][j];
            }
            // Sortowanie iloczynów
            quicksort_float(akumulator, 0, N - 1);
            
            // Użycie algorytmu Kahana do sumowania
            float suma = 0.0;
            float error = 0.0; // Błąd akumulacji
            for (int k = 0; k < N; k++) {
                float y = akumulator[k] - error; // Obliczanie iloczynu z uwzględnieniem błędu
                float t = suma + y;
                error = (t - suma) - y; // Obliczanie nowego błędu
                suma = t; // Aktualizacja sumy
            }
            c[i][j] = suma; // Zapisz wynik
        }
    }
}

void multiply_double_sorted_kahana(double **d, double **e, double **f) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double akumulator[N]; // Tablica do przechowywania iloczynów
            // Mnożenie macierzy
            for (int k = 0; k < N; k++) {
                akumulator[k] = d[i][k] * e[k][j];
            }
            // Sortowanie iloczynów
            quicksort_double(akumulator, 0, N - 1);
            
            // Użycie algorytmu Kahana do sumowania
            double suma = 0.0;
            double error = 0.0; // Błąd akumulacji
            for (int k = 0; k < N; k++) {
                double y = akumulator[k] - error; // Obliczanie iloczynu z uwzględnieniem błędu
                double t = suma + y;
                error = (t - suma) - y; // Obliczanie nowego błędu
                suma = t; // Aktualizacja sumy
            }
            f[i][j] = suma; // Zapisz wynik
        }
    }
}

void multiply_float_sorted_sumator_kahana(float **a, float **b, double **f) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double akumulator[N]; // Tablica do przechowywania iloczynów
            for (int k = 0; k < N; k++) {
                akumulator[k] = a[i][k] * b[k][j];
            }
            // Sortowanie iloczynów
            quicksort_double(akumulator, 0, N - 1);
            
            // Użycie algorytmu Kahana do sumowania
            double suma = 0.0;
            double error = 0.0; // Błąd akumulacji
            for (int k = 0; k < N; k++) {
                double y = akumulator[k] - error; // Obliczanie iloczynu z uwzględnieniem błędu
                double t = suma + y;
                error = (t - suma) - y; // Obliczanie nowego błędu
                suma = t; // Aktualizacja sumy
            }
            f[i][j] = suma; // Zapisz wynik
        }
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

void calculate_error(long double **reference, float **result, 
                     float *min_error_float, float *max_error_float, float *mean_error_float) {
    *max_error_float = 0.0;
    *min_error_float = FLT_MAX;// Ustawienie na maksymalną wartość dla double
    long double total_error = 0.0;
    int count = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Sprawdzenie, czy wartość referencyjna jest różna od zera, aby uniknąć dzielenia przez zero
            if (fabsl(reference[i][j]) > 0) {
                float relative_error = fabsf(result[i][j] - reference[i][j]) / fabsl(reference[i][j]);
                total_error += relative_error;

                // Aktualizacja maksymalnego błędu względnego
                if (relative_error > *max_error_float) {
                    *max_error_float = relative_error;
                }

                // Aktualizacja minimalnego błędu względnego
                if (relative_error < *min_error_float) {
                    *min_error_float = relative_error;
                }

                count++;
            }
        }
    }

    // Obliczanie średniego błędu względnego
    if (count > 0) {
        *mean_error_float = total_error / count;
    } else {
        *mean_error_float = 0.0; // Brak wartości do obliczenia średniego błędu
        *min_error_float = 0.0; // Ustawienie na 0, gdy brak wartości
        *max_error_float = 0.0; // Ustawienie na 0, gdy brak wartości
    }
}

void calculate_error_double(long double **reference, double **result, 
                     double *min_error_double, double *max_error_double, double *mean_error_double) {
    *max_error_double = 0.0;
    *min_error_double = DBL_MAX;// Ustawienie na maksymalną wartość dla double
    long double total_error = 0.0;
    int count = 0;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Sprawdzenie, czy wartość referencyjna jest różna od zera, aby uniknąć dzielenia przez zero
            if (fabsl(reference[i][j]) > 0) {
                long double relative_error = fabsl(result[i][j] - reference[i][j]) / fabsl(reference[i][j]);
                total_error += relative_error;

                // Aktualizacja maksymalnego błędu względnego
                if (relative_error > *max_error_double) {
                    *max_error_double = relative_error;
                }

                // Aktualizacja minimalnego błędu względnego
                if (relative_error < *min_error_double) {
                    *min_error_double = relative_error;
                }

                count++;
            }
        }
    }

    // Obliczanie średniego błędu względnego
    if (count > 0) {
        *mean_error_double = total_error / count;
    } else {
        *mean_error_double = 0.0; // Brak wartości do obliczenia średniego błędu
        *min_error_double = 0.0; // Ustawienie na 0, gdy brak wartości
        *max_error_double = 0.0; // Ustawienie na 0, gdy brak wartości
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

void print_float(float matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void print_double(double matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
}

void print_long_double(long double matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%Lf ", matrix[i][j]);
        }
        printf("\n");
    }
}

void reset(float **c, double **f, long double **g) {      // Reset macierze wynikowe
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c[i][j] = 0;
            f[i][j] = 0;
            g[i][j] = 0;
        }
    }
}

int main() {
    double B = 1000000;
    float MIN_f = -B, MAX_f = B;
    double MIN_d = MIN_f, MAX_d = MAX_f;
    unsigned long long start, end;
    unsigned int lo, hi;
    srand(time(NULL));
    float min_error_float, max_error_float, mean_error_float;
    double min_error_double, max_error_double, mean_error_double;
    long double min_error_long, max_error_long, mean_error_long;
    
    float **a = malloc(N * sizeof(float*));
    float **b = malloc(N * sizeof(float*));
    float **c = malloc(N * sizeof(float*));
    
    double **d = malloc(N * sizeof(double*));
    double **e = malloc(N * sizeof(double*));
    double **f = malloc(N * sizeof(double*));
    
    long double **g = malloc(N * sizeof(long double*));
    long double **reference = malloc(N * sizeof(long double*));
    
    for (int i = 0; i < N; i++) {
        a[i] = malloc(N * sizeof(float));
        b[i] = malloc(N * sizeof(float));
        c[i] = malloc(N * sizeof(float));
        d[i] = malloc(N * sizeof(double));
        e[i] = malloc(N * sizeof(double));
        f[i] = malloc(N * sizeof(double));
        g[i] = malloc(N * sizeof(long double));
        reference[i] = malloc(N * sizeof(long double));
    }
    
    int rounding_mode = fegetround();

    switch (rounding_mode) {
        case FE_TONEAREST:
            printf("Domyślny tryb zaokrąglania: zaokrąglanie do najbliższej liczby\n\n");
            break;
        case FE_DOWNWARD:
            printf("Domyślny tryb zaokrąglania: zaokrąglanie w dół\n\n");
            break;
        case FE_UPWARD:
            printf("Domyślny tryb zaokrąglania: zaokrąglanie w górę\n\n");
            break;
        case FE_TOWARDZERO:
            printf("Domyślny tryb zaokrąglania: zaokrąglanie do zera\n\n");
            break;
        default:
            printf("Nieznany tryb zaokrąglania\n\n");
            break;
    }
    
    // Inicjalizacja macierzy a i b oraz kopiowanie ich do d i e
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            a[i][j] = MIN_f + (MAX_f - MIN_f) * ((float)rand())/RAND_MAX;
            d[i][j] = (double)a[i][j];
            b[i][j] = MIN_f + (MAX_f - MIN_f) * ((float)rand())/RAND_MAX;
            e[i][j] = (double)b[i][j];
            c[i][j] = 0; // Inicjalizacja macierzy c
            f[i][j] = 0; // Inicjalizacja macierzy f
            g[i][j] = 0; // Inicjalizacja macierzy g
        }
    }
    
    /*for (int i = 0; i < N; i++) { // Inicjalizacja macierzy d i e, w celu zbadania wpływu zakresu wartości liczb zmiennoprzecinkowych na błąd, z użyciem najlepszego wariantu
        for (int j = 0; j < N; j++) {
            d[i][j] = MIN_d + (MAX_d - MIN_d) * ((double)rand())/RAND_MAX;
            e[i][j] = MIN_d + (MAX_d - MIN_d) * ((double)rand())/RAND_MAX;
            c[i][j] = 0;
            f[i][j] = 0;
            g[i][j] = 0;
        }
    }*/
    
    multiply_double_sorted_sumator_kahana(d, e, reference);
    
    // Pomiar czasu dla mnożenia naiwnym
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_float(a, b, c);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia float: %llu cykli\n", end - start);
    calculate_error(reference, c, &min_error_float, &max_error_float, &mean_error_float);
    printf("Minimalny błąd względny: %e\n", min_error_float);
    printf("Maksymalny błąd względny: %e\n", max_error_float);
    printf("Średni błąd względny: %e\n\n", mean_error_float);
    
    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_double(d, e, f);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia double: %llu cykli\n", end - start);
    calculate_error_double(reference, f, &min_error_double, &max_error_double, &mean_error_double);
    printf("Minimalny błąd względny: %e\n", min_error_double);
    printf("Maksymalny błąd względny: %e\n", max_error_double);
    printf("Średni błąd względny: %e\n\n", mean_error_double);
    
    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_float_sorted(a, b, c);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia float_sorted: %llu cykli\n", end - start);
    calculate_error(reference, c, &min_error_float, &max_error_float, &mean_error_float);
    printf("Minimalny błąd względny: %e\n", min_error_float);
    printf("Maksymalny błąd względny: %e\n", max_error_float);
    printf("Średni błąd względny: %e\n\n", mean_error_float);
    
    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_double_sorted(d, e, f);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia double_sorted: %llu cykli\n", end - start);
    calculate_error_double(reference, f, &min_error_double, &max_error_double, &mean_error_double);
    printf("Minimalny błąd względny: %e\n", min_error_double);
    printf("Maksymalny błąd względny: %e\n", max_error_double);
    printf("Średni błąd względny: %e\n\n", mean_error_double);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_float_sumator(a, b, f);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia float_sumator: %llu cykli\n", end - start);
    calculate_error_double(reference, f, &min_error_double, &max_error_double, &mean_error_double);
    printf("Minimalny błąd względny: %e\n", min_error_double);
    printf("Maksymalny błąd względny: %e\n", max_error_double);
    printf("Średni błąd względny: %e\n\n", mean_error_double);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_double_sumator(d, e, g);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia double_sumator: %llu cykli\n", end - start);
    calculate_error_long(reference, g, &min_error_long, &max_error_long, &mean_error_long);
    printf("Minimalny błąd względny: %Le\n", min_error_long);
    printf("Maksymalny błąd względny: %Le\n", max_error_long);
    printf("Średni błąd względny: %Le\n\n", mean_error_long);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_float_kahana(a, b, c);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia float_kahana: %llu cykli\n", end - start);
    calculate_error(reference, c, &min_error_float, &max_error_float, &mean_error_float);
    printf("Minimalny błąd względny: %e\n", min_error_float);
    printf("Maksymalny błąd względny: %e\n", max_error_float);
    printf("Średni błąd względny: %e\n\n", mean_error_float);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_double_kahana(d, e, f);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia double_kahana: %llu cykli\n", end - start);
    calculate_error_double(reference, f, &min_error_double, &max_error_double, &mean_error_double);
    printf("Minimalny błąd względny: %e\n", min_error_double);
    printf("Maksymalny błąd względny: %e\n", max_error_double);
    printf("Średni błąd względny: %e\n\n", mean_error_double);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_float_sorted_sumator(a, b, f);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia float_sorted_sumator: %llu cykli\n", end - start);
    calculate_error_double(reference, f, &min_error_double, &max_error_double, &mean_error_double);
    printf("Minimalny błąd względny: %e\n", min_error_double);
    printf("Maksymalny błąd względny: %e\n", max_error_double);
    printf("Średni błąd względny: %e\n\n", mean_error_double);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_double_sorted_sumator(d, e, g);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia double_sorted_sumator: %llu cykli\n", end - start);
    calculate_error_long(reference, g, &min_error_long, &max_error_long, &mean_error_long);
    printf("Minimalny błąd względny: %Le\n", min_error_long);
    printf("Maksymalny błąd względny: %Le\n", max_error_long);
    printf("Średni błąd względny: %Le\n\n", mean_error_long);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_float_sumator_kahana(a, b, f);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia float_sumator_kahana: %llu cykli\n", end - start);
    calculate_error_double(reference, f, &min_error_double, &max_error_double, &mean_error_double);
    printf("Minimalny błąd względny: %e\n", min_error_double);
    printf("Maksymalny błąd względny: %e\n", max_error_double);
    printf("Średni błąd względny: %e\n\n", mean_error_double);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_double_sumator_kahana(d, e, g);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia double_sumator_kahana: %llu cykli\n", end - start);
    calculate_error_long(reference, g, &min_error_long, &max_error_long, &mean_error_long);
    printf("Minimalny błąd względny: %Le\n", min_error_long);
    printf("Maksymalny błąd względny: %Le\n", max_error_long);
    printf("Średni błąd względny: %Le\n\n", mean_error_long);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_float_sorted_kahana(a, b, c);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia float_sorted_kahana: %llu cykli\n", end - start);
    calculate_error(reference, c, &min_error_float, &max_error_float, &mean_error_float);
    printf("Minimalny błąd względny: %e\n", min_error_float);
    printf("Maksymalny błąd względny: %e\n", max_error_float);
    printf("Średni błąd względny: %e\n\n", mean_error_float);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_double_sorted_kahana(d, e, f);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia double_sorted_kahana: %llu cykli\n", end - start);
    calculate_error_double(reference, f, &min_error_double, &max_error_double, &mean_error_double);
    printf("Minimalny błąd względny: %e\n", min_error_double);
    printf("Maksymalny błąd względny: %e\n", max_error_double);
    printf("Średni błąd względny: %e\n\n", mean_error_double);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_float_sorted_sumator_kahana(a, b, f);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia float_sorted_sumator_kahana: %llu cykli\n", end - start);
    calculate_error_double(reference, f, &min_error_double, &max_error_double, &mean_error_double);
    printf("Minimalny błąd względny: %e\n", min_error_double);
    printf("Maksymalny błąd względny: %e\n", max_error_double);
    printf("Średni błąd względny: %e\n\n", mean_error_double);

    reset(c, f, g);
    
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    start = ((unsigned long long)hi << 32) | lo;
    multiply_double_sorted_sumator_kahana(d, e, g);
    __asm__ __volatile__ ("rdtscp" : "=a"(lo), "=d"(hi) : : "%rcx");
    __asm__ __volatile__ ("cpuid" : : : "%rax", "%rbx", "%rcx", "%rdx");
    end = ((unsigned long long)hi << 32) | lo;
    printf("Czas mnożenia double_sorted_sumator_kahana: %llu cykli\n", end - start);
    calculate_error_long(reference, g, &min_error_long, &max_error_long, &mean_error_long);
    printf("Minimalny błąd względny: %Le\n", min_error_long);
    printf("Maksymalny błąd względny: %Le\n", max_error_long);
    printf("Średni błąd względny: %Le\n\n", mean_error_long);

    for (int i = 0; i < N; i++) {
        free(a[i]); free(b[i]); free(c[i]);
        free(d[i]); free(e[i]); free(f[i]);
        free(g[i]); free(reference[i]);
    }
    free(a); free(b); free(c);
    free(d); free(e); free(f);
    free(g); free(reference);

    return 0;
} 
