CC = gcc
MPICC = mpicc
CXX = gxx

CFLAGS = -Wall -g -fopenmp 
CFLAGS = -fopenmp -O3 
CXXFLAGS = -Wall -g -fopenmp

OBJS = mul_matrix_from_file.o matrix.o 

.SUFFIXES: .o .c .cpp .cxx

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<


.PHONY: all clean
#all:    multiply_matrix
all: mul_matrix
#all: mul_matrix_from_file


clean:
	$(RM) *.o

mul_matrix_from_file: mul_matrix_from_file.c
#	$(CC) $(CFLAGS) -o $@ $^ $(BLASLIB)
	$(MPICC) $(CFLAGS) -o multiply_matrix $^

multiply_matrix: multiply_matrix.cpp
	$(MPICC) $(CFLAGS) -o multiply_matrix $~
	
mul_matrix: $(OBJS)
	$(MPICC) $(CFLAGS) -o mul_matrix $(OBJS)