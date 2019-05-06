CFLAGS = -fPIC -Wall -Wno-unused-result -O3 -fopenmp
LFLAGS = -lm -fopenmp -lgsl -lgslcblas
CC = mpicc
OBJS = base.o cspline.o md_ghost.o forces.o potential.o
SHAREDLIB = libfmd.so
all: shared
shared: $(SHAREDLIB)
$(SHAREDLIB): $(OBJS)
	$(CC) $^ -shared -o $@ $(LFLAGS)
base.o: base.h md_ghost.h forces.h
forces.o: forces.h base.h cspline.h md_ghost.h
md_ghost.o: base.h md_ghost.h
potential.o: base.h cspline.h
cspline.o: cspline.h
clean:
	rm -f *.o libfmd.so
