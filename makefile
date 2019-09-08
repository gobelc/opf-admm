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
	rm solve_opf.o solve_opf ./0/*.dat ./1/*.dat ./2/*.dat ./3/*.dat ./4/*.dat ./5/*.dat ./6/*.dat ./7/*.dat ./0/*.txt ./1/*.txt ./2/*.txt ./3/*.txt ./4/*.txt ./5/*.txt ./6/*.txt ./7/*.txt

.PHONY : run
run :
	mpirun -np 8 ./solve_opf