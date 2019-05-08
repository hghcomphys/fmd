CFLAGS = -fPIC -Wall -Wno-unused-result -O3 -fopenmp
LFLAGS = -lm -fopenmp -lgsl -lgslcblas
CC = mpicc
OBJS = base.o cspline.o md_ghost.o forces.o potential.o
SHAREDLIB = libfmd.so
BASE_H = base.h potential.h

.PHONY: all clean
all: shared
shared: $(SHAREDLIB)
$(SHAREDLIB): $(OBJS)
	$(CC) $^ -shared -o $@ $(LFLAGS)
base.o: $(BASE_H) md_ghost.h forces.h
forces.o: forces.h $(BASE_H) cspline.h md_ghost.h
md_ghost.o: $(BASE_H) md_ghost.h
potential.o: $(BASE_H) cspline.h potential.h
cspline.o: cspline.h
clean:
	rm -f *.o libfmd.so
