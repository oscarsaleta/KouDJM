CFLAGS= -std=c11 -Wall -O3
LFLAGS= -lm -lgmp -lmpfr -lmpc

target: all

all: f1 f2 f3 resultant

# Main executable
f1: f1.o koujdm_lib.o
	gcc $(CFLAGS) f1.o koujdm_lib.o -o f1 $(LFLAGS)

f2: f2.o koujdm_lib.o
	gcc $(CFLAGS) f2.o koujdm_lib.o -o f2 $(LFLAGS)

f3: f3.o koujdm_lib.o
	gcc $(CFLAGS) f3.o koujdm_lib.o -o f3 $(LFLAGS)

resultant: resultant.o koujdm_lib.o
	gcc $(CFLAGS) resultant.o koujdm_lib.o -o resultant $(LFLAGS)

# Objects compilation
f1.o: f1.c
	gcc $(CFLAGS) -c f1.c $(LFLAGS)

f2.o: f2.c
	gcc $(CFLAGS) -c f2.c $(LFLAGS)

f3.o: f3.c
	gcc $(CFLAGS) -c f3.c $(LFLAGS)

resultant.o: resultant.c
	gcc $(CFLAGS) -c resultant.c $(LFLAGS)

koujdm_lib.o: koujdm_lib.c koujdm_lib.h
	gcc $(CFLAGS) -c koujdm_lib.c $(LFLAGS)


# Run configuration
runf1: f1
	./f1 0.1 0.2 3 50 33.3333333333333333 0.5 1 0.3 14 12 4

runf2: f2
	./f2 0.1 0.2 3 50 33.3333333333333333 0.5 1 0.2 0.3 14 12 4

runf3: f3
	./f3 500 0.1 0.2 3 50 100/3 0.5 1 0.3 10 2

runresultant: resultant
	./resultant 0.1 0.2 3 50 33.33333333333 0.5


# Cleaning directives
clean:
	rm *o f1 f2 f3 resultant
