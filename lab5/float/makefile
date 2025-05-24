all: fpu

fpu: fpu.c
	gcc -O2 -mfpmath=387 -march=native -o fpu fpu.c -lm

run: fpu
	./fpu

asm: fpu.c
	gcc -O2 -S -o fpu.s fpu.c

disasm: fpu
	objdump -d fpu > fpu_disasm.txt

clean:
	rm -f fpu fpu.s fpu_disasm.txt

.PHONY: clean run asm disasm
