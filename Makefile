CC = g++
LINK = -o
CXX_VERSION = -std=c++14 -pthread
OPTIMIZATION = -O3 -c
OPTIMIZATION_LINK = -O3
CXXFLAG = -Wall $(CXX_VERSION)
SRC = ./src
LIB = ./lib

OBJS = main.o solver.o hungarian.o history.o load_test.o
PROG = sop_solver

$(PROG): $(OBJS)
	$(CC) $(CXXFLAG) $(OPTIMIZATION_LINK) $(LINK) $(PROG) $^ -lgmp -lgmpxx

main.o: $(SRC)/main.cpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $<

solver.o: $(LIB)/solver.cpp $(LIB)/hash.hpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $<

hungarian.o: $(LIB)/hungarian.cpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $<

history.o: $(LIB)/history.cpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $<

load_test.o: $(LIB)/load_test.cpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $<

clean:
	rm -rf *.o sop_solver