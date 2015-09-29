CXX = g++
CXXFLAGS = -g -Wall -O
LDADD =

CPX_INCDIR  = /opt/CPLEX_Studio/cplex/cplex/include
CPX_LDPATH  = /opt/CPLEX_Studio/cplex/cplex/lib/x86-64_linux/static_pic
CPX_LDADD = -lcplex -lm -pthread

OBJ = main.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(CPX_INCDIR) -c $^ -o $@

main: $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o main -L$(CPX_LDPATH) $(CPX_LDADD)

clean:
	rm -rf $(OBJ) main

.PHONY: clean