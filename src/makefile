PREFIX = /usr
CC = mpicc
SHAREDLIB = libfmd.so
CFLAGS = -fPIC -Wall -Wno-unused-result -O3 -fopenmp
LFLAGS = -fopenmp -lgsl -lgslcblas

OBJS ::= $(patsubst %.c,%.o,$(wildcard *.c))

.PHONY: all shared clean install uninstall

all: shared

shared: $(SHAREDLIB)

$(SHAREDLIB): $(OBJS)
	$(CC) $(OBJS) -shared -o $@ $(LFLAGS)

##### the following lines can be updated by using "gcc -MM *.c"
array.o: array.c array.h config.h
base.o: base.c base.h config.h potential.h types.h md_ghost.h forces.h
cspline.o: cspline.c cspline.h config.h
eam.o: eam.c eam.h config.h forces.h types.h cspline.h base.h potential.h \
 list.h
forces.o: forces.c forces.h config.h eam.h types.h cspline.h lj.h morse.h \
 base.h potential.h md_ghost.h list.h
list.o: list.c list.h config.h
lj.o: lj.c lj.h config.h potential.h types.h base.h list.h forces.h
md_ghost.o: md_ghost.c base.h config.h potential.h types.h md_ghost.h
morse.o: morse.c morse.h config.h base.h potential.h types.h forces.h \
 list.h
potential.o: potential.c potential.h config.h types.h base.h array.h \
 list.h eam.h forces.h cspline.h
timer.o: timer.c timer.h config.h base.h potential.h types.h
#####

clean:
	rm -f $(OBJS) $(SHAREDLIB)

install: $(SHAREDLIB)
	mkdir -p $(DESTDIR)$(PREFIX)/lib
	mkdir -p $(DESTDIR)$(PREFIX)/include
	cp fmd.h $(DESTDIR)$(PREFIX)/include/
	cp $(SHAREDLIB) $(DESTDIR)$(PREFIX)/lib/

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/include/fmd.h
	rm -f $(DESTDIR)$(PREFIX)/lib/$(SHAREDLIB)
