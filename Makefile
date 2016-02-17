CFLAGS=-O3 -std=gnu11 -fPIC
LDFLAGS=-shared
CC=clang

SOURCES=engine/main.c engine/jacobi.c engine/sor.c
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
	-rm *.o solver.so
