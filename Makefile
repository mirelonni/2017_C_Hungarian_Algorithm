CC = gcc
CFLAGS = -Wall -g

build:solve.h main.c
	$(CC) main.c -o main $(CFLAGS)

clean:
	rm -f main
