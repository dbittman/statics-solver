CFLAGS=-Ofast -std=gnu11 -fPIC -march=native -Wall -Wextra -msse3 -ffast-math -g -mavx
LDFLAGS=-shared -lpthread
CC=gcc

SOURCES=engine/main.c engine/sor.c engine/jacobi.c
OBJECTS=$(SOURCES:.c=.o)
DEPS=$(SOURCES:.c=.d)

all: engine/solver.so solve

.PHONY: all

solve: frontend/solve
	@cp frontend/solve solve
	@chmod a+x solve

engine/solver.so: $(OBJECTS)
	clang $(LDFLAGS) -Wl,-soname,engine/solver.so -o engine/solver.so $(OBJECTS) -lm

-include $(DEPS)

%.o: %.c
	$(CC) $(CFLAGS) -MF $*.d -MMD -c $< -o $@

clean:
	-rm $(OBJECTS) engine/solver.so solve
