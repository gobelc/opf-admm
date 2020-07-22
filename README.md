# opf-admm [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3956696.svg)](https://doi.org/10.5281/zenodo.3956696)

Distributed Optimal Power Flow with ADMM via Message Passing Interface (MPI) implemented in C++.

[Report (Spanish)](https://iie.fing.edu.uy/~gbelcredi/hpc/informe.pdf)

![alt text](https://github.com/gobelc/opf-admm/blob/master/residuo.png)

Based on Q. Peng and S. H. Low, "Distributed Optimal Power Flow Algorithm for Radial Networks, I: Balanced Single Phase Case," in IEEE Transactions on Smart Grid, vol. 9, no. 1, pp. 111-121, Jan. 2018.

# Requirements
- Matlab
- CVX
- mpich

# Setup
For the setup you need to:
- edit the path to the opf-admm directory in the file NetworkModel.hpp
- edit the path to the matlab binaries, i.e "/usr/local/MATLAB/R2017a/bin/matlab" in files m/update_x.sh and m/update_y.sh
- for each node you can modify the SOCP optimization problem, editing the files 0,1...,7/update_x.m and 0,1...,7/update_y.m

The default program works with 8 buses and 8 MPI processors, this can be configured changing the file NetworkModel.hpp and adjusting the makefile.

## Compiling
make

## Run
make run

## Clean
You can delete previous data with:

make clean

## Example
8-bus configuration example:
![alt text](https://github.com/gobelc/opf-admm/blob/master/grafo.png)

Active power injected:
![alt text](https://github.com/gobelc/opf-admm/blob/master/potencia_inyectada.png)
