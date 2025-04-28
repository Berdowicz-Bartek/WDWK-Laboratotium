all: fpu

fpu: fpu.c
	gcc -O0 -ffloat-store -mno-sse -o fpu fpu.c -lm -msse

run: fpu
	./fpu

asm: fpu.c
	gcc -O0 -S -o fpu.s fpu.c

disasm: fpu
	objdump -d fpu > fpu_disasm.txt

clean:
	rm -f fpu fpu.s fpu_disasm.txt

.PHONY: clean run asm disasm
