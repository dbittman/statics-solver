CFLAGS=-O3 -std=gnu11 -fPIC -flto
LDFLAGS=-shared -flto
CC=gcc

solver.so: main.o jacobi.o
	gcc $(LDFLAGS) -Wl,-soname,solver.so -o solver.so jacobi.o main.o -lm

main.o: main.c grid.h

jacobi.o: jacobi.c grid.h

clean:
	-rm *.o solver.so
