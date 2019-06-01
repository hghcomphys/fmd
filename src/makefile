CC = mpicc
SHAREDLIB = libfmd.so
CFLAGS = -fPIC -Wall -Wno-unused-result -O3 -fopenmp
LFLAGS = -lm -fopenmp -lgsl -lgslcblas

OBJS ::= $(patsubst %.c,%.o,$(wildcard *.c))

.PHONY: all shared clean
all: shared
shared: $(SHAREDLIB)
$(SHAREDLIB): $(OBJS)
	$(CC) $(OBJS) -shared -o $@ $(LFLAGS)

array.o: array.c array.h
base.o: base.c base.h potential.h md_ghost.h forces.h
cspline.o: cspline.c cspline.h
forces.o: forces.c forces.h base.h potential.h cspline.h md_ghost.h \
 list.h
list.o: list.c list.h
md_ghost.o: md_ghost.c base.h potential.h md_ghost.h
potential.o: potential.c potential.h base.h cspline.h array.h list.h

clean:
	rm -f *.o libfmd.so
