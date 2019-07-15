# Specify compiler
CC=gcc

# Specify linker
LINK=g++

.PHONY : all
all : solve_opf

# Link the object files into a binary
solve_opf : solve_opf.o
	$(LINK) -o solve_opf solve_opf.o

# Compile the source files into object files
solve_opf.o : solve_opf.cpp
	$(CC) -c -std=c++11 solve_opf.cpp -lm -o solve_opf.o 

# Clean target
.PHONY : clean
clean :
	rm solve_opf.o solve_opf 

.PHONY : run
run :
	./solve_opf 