CC = g++
LINK = $(CC) -o
CXX_VERSION = -std=c++11
OPTIMIZATION = -O3 -Og -pg
CXXFLAG = -c -g $(CXX_VERSION) -Wall
SRC = ./src
LIB = ./lib

OBJS = main.o solver.o hungarian.o history.o
PROG = sop_solver

$(PROG): $(OBJS)
	$(LINK) $(PROG) $^

main.o: $(SRC)/main.cpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $<

solver.o: $(LIB)/solver.cpp $(LIB)/hash.hpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $<

hungarian.o: $(LIB)/hungarian.cpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $<

history.o: $(LIB)/history.cpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $<

clean:
	rm -rf *.o sop_solver