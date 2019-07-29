MPI_COMPILE_FLAGS = $(shell mpic++ --showme:compile)
MPI_LINK_FLAGS = $(shell mpic++ --showme:link)

# Specify compiler
CC=mpic++

.PHONY : all
all : solve_opf

# Compile the source files into object files
solve_opf : solve_opf.cpp
	$(CC) $(MPI_COMPILE_FLAGS) -std=c++11 solve_opf.cpp $(MPI_LINK_FLAGS) -o solve_opf

# Clean target
.PHONY : clean
clean :
	rm solve_opf.o solve_opf 

.PHONY : run
run :
	mpirun -np 4 ./solve_opf 