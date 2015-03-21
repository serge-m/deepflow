CC=gcc

CFLAGS=-Wall -g -O3
LDFLAGS=-g -Wall -O3
LIBFLAGS=-lm -ljpeg -lpng
LIBAFLAGS=-static /usr/lib64/libjpeg.a /usr/lib64/libpng.a /usr/lib64/libz.a /usr/lib64/libm.a

all: deepflow

deepflow: deepflow.o image.o io.o opticalflow_aux.o opticalflow.o solver.o
	$(CC) $(LDFLAGS) $(LIBFLAGS) -o $@ $^

deepflow-static: deepflow.o image.o io.o opticalflow_aux.o opticalflow.o solver.o
	$(CC) $(LIBFLAGS) -o $@ $^ $(LIBAFLAGS)

%.o: %.c
	$(CC) -o $@ $(CFLAGS) -c $+ 

clean:
	rm -f *.o deepflow
