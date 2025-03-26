main: quicksort.o main.o
	gcc -no-pie main.o quicksort.o -o main

test.o: quicksort.c
	gcc -Wall -Wextra -c quicksort.c -o quicksort.o

asm.o: main.s
	as -o main.o main.s
	
run: main
	./main

clean:
	rm -f quicksort.o main.o main
	
.PHONY: clean run
