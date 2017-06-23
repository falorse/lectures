CC = gcc
MPICC = mpicc
CXX = gxx

//CFLAGS = -Wall -g -fopenmp 
CFLAGS = -fopenmp -O3 
CXXFLAGS = -Wall -g -fopenmp

OBJS = matrix.o

.SUFFIXES: .o .c .cpp .cxx

.c.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<


.PHONY: all clean
all:    mul_matrix_from_file

clean:
	$(RM) *.o

mul_matrix_from_file: mul_matrix_from_file.c
#	$(CC) $(CFLAGS) -o $@ $^ $(BLASLIB)
	$(MPICC) $(CFLAGS) -o multiply_matrix $^

