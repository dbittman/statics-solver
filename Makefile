CFLAGS=-Ofast -std=gnu11 -flto -fPIC
LDFLAGS='-lm' -flto -shared
CC=gcc

solver.so: main.o jacobi.o
	gcc $(LDFLAGS) -o solver.so main.o jacobi.o

main.o: main.c grid.h

jacobi.o: jacobi.c grid.h

clean:
	rm *.o solver.sm
