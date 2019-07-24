A = csvread("/home/gonzalo/workspace/opf-admm/A.csv");
x = csvread("/home/gonzalo/workspace/opf-admm/x.csv");
y = csvread("/home/gonzalo/workspace/opf-admm/y.csv");
mu = csvread("/home/gonzalo/workspace/opf-admm/mu.csv");


%%% Init variables OPF
n = size(y);
rho=0.25;

% Problem OPF: 
cvx_begin
    variable y(n,1)
    minimize(-mu'*y  + .5*rho*sum_square(x-y));
    A*y == 0;
cvx_end

mu = mu + rho*(x - y);
csvwrite('/home/gonzalo/workspace/opf-admm/y.csv',y)
csvwrite('/home/gonzalo/workspace/opf-admm/mu.csv',mu)

