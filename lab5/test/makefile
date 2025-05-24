all: test simd

test: test.c
	gcc -O2 -mfpmath=387 -march=native -mno-sse -mno-sse2 -mno-sse3 -mno-ssse3 -mno-sse4.1 -mno-sse4.2 -mno-avx -mno-avx2 -o test test.c

simd: test.c
	gcc -O2 -mavx2 -mfma -march=native -Wall -Wextra -o test_simd test.c

run: test
	./test

runSIMD: simd
	./test_simd
	
asm: test.c
	gcc -O2 -S -o test.s test.c

disasm: test
	objdump -d test > test_disasm.txt
	
asmSIMD: test.c
	gcc -O2 -mavx2 -mfma -march=native -Wall -Wextra -S -o test_simd.s test.c

disasmSIMD: test_simd
	objdump -d test_simd > test_simd_disasm.txt

clean:
	rm -f test test.s test_disasm.txt test_simd test_simd.s test_simd_disasm.txt

.PHONY: clean run
