#PATH to where RIsearch tables are to be found
MATPATH = \"/home/projects/users/anthon/projects/2024_risearch/temp_par/risearch/tables\"
CC = gcc-15
CFLAGS += -fopenmp -Wall -I../libdivsufsort-2.0.1/build/include -DRISEARCH2 -DMATPATH=$(MATPATH) -DRISVERSION=2
CLIBS += -lz -lpcre -lm

all: ../bin/risearch2.x ../bin/risearch2.dbg.x

SOURCES := $(filter-out dump_dsm.c,$(wildcard *.c))

../bin/risearch2.x: $(SOURCES) *.h
	$(CC) $(CFLAGS) -O3 $(SOURCES) ../libdivsufsort-2.0.1/build/lib/libdivsufsort64.a $(CLIBS) -o $@

../bin/risearch2.dbg.x: $(SOURCES) *.h
	$(CC) $(CFLAGS) -ggdb -g3 -O0 -DDEBUG -DVERBOSE=2 $(SOURCES) ../libdivsufsort-2.0.1/build/lib/libdivsufsort64.a $(CLIBS) -o $@

# Profiling build: optimized with debug symbols for flamegraph
../bin/risearch2.prof.x: $(SOURCES) *.h
	$(CC) $(CFLAGS) -g -O3 $(SOURCES) ../libdivsufsort-2.0.1/build/lib/libdivsufsort64.a $(CLIBS) -o $@

.PHONY: clean

clean:
	rm -f ../bin/risearch2.x ../bin/risearch2.dbg.x
