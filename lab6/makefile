OUT ?= results.csv
PERF ?= perf.data
PERF_IN ?= cache_gcc_O0

gcc_O0: cache.c
	gcc -O0 -Wall -o cache_gcc_O0 cache.c

gcc_O2: cache.c
	gcc -O2 -Wall -o cache_gcc_O2 cache.c

clang_O0: cache.c
	clang -O0 -Wall -o cache_clang_O0 cache.c

clang_O2: cache.c
	clang -O2 -Wall -o cache_clang_O2 cache.c

all: gcc_O0 gcc_O2 clang_O0 clang_O2

run_gcc_O0:
	./cache_gcc_O0 > $(OUT)

run_gcc_O2:
	./cache_gcc_O2 > $(OUT)

run_clang_O0:
	./cache_clang_O0 > $(OUT)

run_clang_O2:
	./cache_clang_O2 > $(OUT)

perf:
	perf record -o $(PERF) -- ./$(PERF_IN)

clean:
	rm -f cache_gcc_O0 cache_gcc_O2 cache_clang_O0 cache_clang_O2 *.csv perf.data perf.data.old

.PHONY: all gcc_O0 gcc_O2 clang_O0 clang_O2 run perf clean

