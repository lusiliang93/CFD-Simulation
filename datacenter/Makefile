CC=gcc
CFLAGS = -Wall -pedantic
.SUFFIXES: .o .c
.c.o: ; $(CC) -c $(CFLAGS) $*.c
OBJ = init.o boundary.o uvp.o main.o

run: $(OBJ)
	$(CC) $(CFLAGS) -o run $(OBJ) -lm

init.o :init.h
boundary.o: boundary.h
uvp.o: uvp.h
main.o: init.h boundary.h uvp.h
