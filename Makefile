CC         = gcc
CFLAGS     = -c -lm -Wall
GCFLAGS    = ./Reorderings/hsl_mc73/hsl_mc73d.o ./Reorderings/hsl_mc73/libhsl_mc73.a -lblas  -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lm -lgfortran
SOURCES    = ./CommonFiles/blas.c          \
	     ./CommonFiles/matrix.c        \
	     ./CommonFiles/linked_list.c   \
	     ./CommonFiles/graph.c         \
	     ./CommonFiles/graph_parallel.c\
	     ./Reorderings/rcm.c           \
	     ./Reorderings/rcm_unordered.c \
             ./Reorderings/sloan.c         \
             ./Reorderings/spectral.c      \
	     program.c
OBJECTS    = $(SOURCES:rail_5177.mtx.c=.o)
EXECUTABLE = program 

all: 	$(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(GCFLAGS) -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

run:
	./$(EXECUTABLE) Matrices/rail_5177.mtx


