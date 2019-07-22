A = csvread("/home/olaznog/workspace/opf-admm/A.csv");
x = csvread("/home/olaznog/workspace/opf-admm/x.csv");
y = csvread("/home/olaznog/workspace/opf-admm/y.csv");
mu = csvread("/home/olaznog/workspace/opf-admm/mu.csv");


%%% Init variables OPF
n = size(y);
rho=.1;

% Problem OPF: 
cvx_begin
    variable y(n)
    minimize(-mu'*y  + .5*rho*sum_square(x-y));
    A*y == 0;
cvx_end

csvwrite('/home/olaznog/workspace/opf-admm/y.csv',y)
