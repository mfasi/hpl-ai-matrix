SHELL = /bin/sh
CC=mpicc
CFLAGS=-Wall -Ofast
CFLAGSHPC=-DMKL_ILP64 -xsse4.2 -axCORE-AVX512
CLIBS=-lmkl_scalapack_ilp64 -lmkl_blacs_openmpi_ilp64 \
	-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 \
	-lpthread -lm -lgsl
CMAKE_PATH=-L/opt/intel/compilers_and_libraries_2019/linux/mkl/lib/intel64_lin	\
	-L/opt/intel/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin
CMAKE_INCLUDE=-I/opt/intel/compilers_and_libraries_2019/linux/mkl/include
DEPS=print_util.h matnorm.h matgen.h
OBJ=print_util.o matnorm.o matgen.o

all: test_lowprec test_norminf

test_lowprec:
	gcc -c print_util.c $(CFLAGS)
	gcc -c matnorm.c $(CFLAGS)
	gcc -c test_lowprec.c $(CFLAGS)
	gcc -Wall -Ofast -lm -lgsl print_util.o matnorm.o test_lowprec.o -o test_lowprec

test_norminf: test_norminf.o $(OBJ)
	$(CC) $(CFLAGS) $(CFLAGSHPC) -o $@ $^ $(CLIBS)

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(CFLAGSHPC)

clean:
	rm test_lowprec test_norminf *.o
	rm *.e* *.o* *~
