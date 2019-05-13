CFLAGS = -fPIC -Wall -Wno-unused-result -O3 -fopenmp
LFLAGS = -lm -fopenmp -lgsl -lgslcblas
CC = mpicc
OBJS = base.o cspline.o md_ghost.o forces.o potential.o array.o list.o
SHAREDLIB = libfmd.so
BASE_H = base.h potential.h

.PHONY: all shared clean
all: shared
shared: $(SHAREDLIB)
$(SHAREDLIB): $(OBJS)
	$(CC) $^ -shared -o $@ $(LFLAGS)
base.o: $(BASE_H) md_ghost.h forces.h
forces.o: forces.h $(BASE_H) cspline.h md_ghost.h
md_ghost.o: $(BASE_H) md_ghost.h
potential.o: $(BASE_H) cspline.h potential.h array.h list.h
cspline.o: cspline.h
array.o: array.h
list.o: list.h
clean:
	rm -f *.o libfmd.so
