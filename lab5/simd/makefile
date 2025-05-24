all: simd

simd: simd.c
	gcc -O2 -mavx2 -mfma -march=native -Wall -Wextra -o simd simd.c

run: simd
	./simd

asm: simd.c
	gcc -O2 -mavx2 -mfma -march=native -Wall -Wextra -S -o simd.s simd.c

disasm: simd
	objdump -d simd > simd_disasm.txt

clean:
	rm -f simd simd.s simd_disasm.txt

.PHONY: all run asm disasm clean

