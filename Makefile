CC=g++

CFLAGS=-Wall -g -O3 -msse4
LDFLAGS=-g -Wall -O3 -msse4
LIBFLAGS=-lm -ljpeg -lpng
LIBAFLAGS=-static /usr/lib/libjpeg.a /usr/lib/libpng.a /usr/lib/libz.a /usr/lib/libm.a

all: fastdeepflow

fastdeepflow: deepflow.o image.o io.o opticalflow_aux.o opticalflow.o solver.o
	$(CC) $(LDFLAGS)  -o $@ $^ $(LIBFLAGS)

fastdeepflow-static: deepflow.o image.o io.o opticalflow_aux.o opticalflow.o solver.o
	$(CC) $(LIBFLAGS) -o $@ $^ $(LIBAFLAGS)

%.o: %.c
	$(CC) -o $@ $(CFLAGS) -c $+ 

clean:
	rm -f *.o deepflow
