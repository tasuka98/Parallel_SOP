CC = g++
LINK = $(CC) -o
CXX_VERSION = -std=c++11
CXXFLAG = -c -g $(CXX_VERSION) -Wall -Werror
SRC = ./src
LIB = ./lib

OBJS = main.o solver.o hungarian.o
PROG = sop_solver

$(PROG): $(OBJS)
	$(LINK) $(PROG) $^

main.o: $(SRC)/main.cpp
	$(CC) $(CXXFLAG) $<

solver.o: $(LIB)/solver.cpp
	$(CC) $(CXXFLAG) $<

hungarian.o: $(LIB)/hungarian.cpp
	$(CC) $(CXXFLAG) $<

clean:
	rm -rf *.o sop_solver