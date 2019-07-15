clc
clear all
close all


%%% Init variables OPF

R1 = .2
R2 = .2
R3 = .2
X1 = .01
X2 = .01
X3 = .01

A = [1, -1, 0, 0, 2*R1 2*X1, -(R1^2 + R2^2),0,0,0,0,0,0,0,0;
    0,0,1,0,-1,0,0,1,0,-R1,1,1,0,-R2,0;
    0,0,0,1,0,-1,0,0,1,-X1,0,0,1,-X2,1
]




% Problem OPF: 
disp('Computing optimal solution for 2nd formulation...');
cvx_begin
    variable x2(n)
    variable Y(n,n) diagonal
    minimize( matrix_frac(A*x2 + b , eye(m) + B*Y*B') )
    x2 >= 0;
    Y == diag(x2);
cvx_end
opt2 = cvx_optval;
