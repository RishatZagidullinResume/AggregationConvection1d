GCC = g++

LIB  = -lgomp -lm
COPT = -fopenmp -O3
GCC_COPT = -pthread -std=c++11
VPATH = ../:../solver

SOURCES = solver.cpp main.cpp
OBJECTS = $(SOURCES:.cpp=_cpp.o)

all: solver.exe

solver.exe: $(OBJECTS)
	$(GCC) $^ -o $@ $(GCC_COPT) $(COPT) $(LIB)

%_cpp.o: %.cpp
	$(GCC) -c $< -o $@ $(GCC_COPT) $(COPT) 

clean:
	rm -rf *.o
	rm -rf *.exe
