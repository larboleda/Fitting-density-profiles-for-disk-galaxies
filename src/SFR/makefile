CC=gcc
CFLAGS=-g -Wall -c -I/$HOME/local/include
LDFLAGS=-lm -L/$HOME/local/lib -lgsl -lgslcblas 
SOURCES=tasa.c 
OBJECTS=$(SOURCES:%.c=%.o)

all: main

main: $(OBJECTS) 
	$(CC)  $(OBJECTS) -o $@.x $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o *.x *~
