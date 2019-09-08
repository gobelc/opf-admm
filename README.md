# opf-admm

Distributed Optimal Power Flow with ADMM via message passing

# Setup
For the setup you need to:
- edit the path to the opf-admm directory in the file NetworkModel.hpp
- edit the path to the matlab binaries, i.e "/usr/local/MATLAB/R2017a/bin/matlab" in files m/update_x.sh and m/update_y.sh

The default program works with 8 buses and 8 MPI processors, this can be configured changing the file NetworkModel.hpp and adjusting the makefile.

# Compiling
make

# Run
make run

You can delete previous data with:

make clean
