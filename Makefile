CFLAGS=-Ofast -std=gnu11
LDFLAGS='-lm'
CC=gcc

main: main.o

main.o: main.c

clean:
	rm main.o main
