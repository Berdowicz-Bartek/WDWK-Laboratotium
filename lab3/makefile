all: matrix

matrix: matrix.c
	gcc -O2 -o matrix matrix.c

run: matrix
	./matrix

clean:
	rm -f matrix

.PHONY: clean run
