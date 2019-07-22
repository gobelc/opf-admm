clear all

A = csvread("/home/olaznog/workspace/opf-admm/A.csv");
x_old = csvread("/home/olaznog/workspace/opf-admm/x.csv");
y = csvread("/home/olaznog/workspace/opf-admm/y.csv");
mu = csvread("/home/olaznog/workspace/opf-admm/mu.csv");

%%% Init variables OPF
n = size(y);
c = zeros(n);
c(3)=1;

cp = zeros(n);
cp(5)=200;
rho=.1;

% Problem OPF: 
cvx_begin
    variable x(n)
    minimize(square(c'*x) + mu'*x + .5*rho*sum_square(x-y));
    x(1)==1.;
    x(3)<=200;
    x(3)>=100
    x(2)<=1.1;
    x(2)>=.9;
    %square(x(5)) + square(x(6)) <= x(2)*x(7);
cvx_end
disp(norm(x-y))
csvwrite('/home/olaznog/workspace/opf-admm/x.csv',x)
