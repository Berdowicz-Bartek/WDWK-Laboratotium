all: matrix

matrix: matrix.c
	gcc -O2 -o matrix matrix.c

run: matrix
	./matrix
	
asm: matrix.c
	gcc -O2 -S -o matrix.s matrix.c

disasm: matrix
	objdump -d matrix > matrix_disasm.txt

clean:
	rm -f matrix matrix.s matrix_disasm.txt

.PHONY: clean run asm disasm
