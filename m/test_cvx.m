clc
clear all
close all

A = csvread("../A.dat");
x = csvread("../x.dat");
y = csvread("../y.dat");
mu = csvread("../mu.dat");


%%% Init variables OPF
n = size(y)
B = A*y;

mu = rand(n)
c = ones(n);
c(3)=1;
%cost = mu'*y ;
rho=1;
y = ones(n)
n_iter = 10;
y_vector = zeros(13,n_iter);
x_vector = zeros(13,n_iter);
residual=zeros(n_iter,1);
p_gen=zeros(n_iter,1);
for k=1:n_iter

    % Problem OPF: 
    disp('Computing optimal solution for 2nd formulation...');
    cvx_begin
        variable x(n)
        minimize(c'*x + mu'*x + .5*rho*sum_square(x-y));
        x(3)>=200;
    cvx_end
    fprintf(1,'Done! \n');

    % Problem OPF: 
    disp('Computing optimal solution for 2nd formulation...');
    cvx_begin
        variable y(n)
        minimize(-mu'*y  + .5*rho*sum_square(x-y));
        A*y == 0;
    cvx_end
    fprintf(1,'Done! \n');
    y_vector(:,k) = y;
    x_vector(:,k) = x;
    residual(k) = norm(x-y);
    p_gen(k) = x(3);
    mu = mu +rho*(x - y);
    
    % Display results
%     disp('------------------------------------------------------------------------');
%     disp('The computed optimal solution is: ');
%     disp(x)
%     fprintf('Q gen = %d kVA\n' , x(4))
%     fprintf('P gen = %d kW\n' , x(3))

    % Display results
     disp('------------------------------------------------------------------------');
     disp('The computed optimal solution is: ');
     disp(y)
     fprintf('Q gen = %d VA\n' , y(4))
     fprintf('P gen = %d W\n' , y(3))
end

disp('------------------------------------------------------------------------');
disp('The computed optimal solution is: ');
disp(y)
fprintf('Q gen = %d kVA\n' , y(4))
fprintf('P gen = %d kW\n' , y(3))
