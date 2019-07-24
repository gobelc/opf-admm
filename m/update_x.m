clear all

A = csvread("/home/gonzalo/workspace/opf-admm/A.csv");
x_old = csvread("/home/gonzalo/workspace/opf-admm/x.csv");
y = csvread("/home/gonzalo/workspace/opf-admm/y.csv");
mu = csvread("/home/gonzalo/workspace/opf-admm/mu.csv");

%%% Init variables OPF
n = size(y,1);
c = zeros(n,1);
c(3)=1;

cp = zeros(n);
cp(5)=200;
rho=0.25;

n_childs = floor((n-7)/3);
B = zeros(3,n);
B(1,5)=2;
B(2,6)=2;
B(3,2)=-1;
B(3,7)=1;

D = zeros(n,1); 
E = zeros(n,1);

D(1,1) = 1;
E(7,1) = 1;

% Problem OPF: 
cvx_begin
    variable x(n,1)
    minimize(square(c'*x) + mu'*x + .5*rho*sum_square(x-y));
    x(1)==1.;
    x(3)<=200;
    x(3)>=100
    x(2)<=1.1;
    x(2)>=.9;
    norm(B*x,2) <= x(1) + x(7);
cvx_end
fprintf("Residuo: %f",norm(x-y,2))
csvwrite('/home/gonzalo/workspace/opf-admm/x.csv',x)
