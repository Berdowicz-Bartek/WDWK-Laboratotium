all:
	as multi.s -o multi.o 
	ld multi.o -o multi -lc -dynamic-linker /lib64/ld-linux-x86-64.so.2
	./multi
