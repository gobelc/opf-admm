function [] = update_y(arg1)

%fprintf(strcat(arg1,"A.csv"))

A = csvread(strcat(arg1,"A.csv"));
x = csvread(strcat(arg1,"x.csv"));
y = csvread(strcat(arg1,"y.csv"));
mu = csvread(strcat(arg1,"mu.csv"));


%%% Init variables OPF
n = size(y);
rho=0.1;

% Problem OPF: 
cvx_begin quiet
    variable y(n,1)
    minimize(-mu'*y  + .5*rho*power(2,norm(x-y,2)));
    A*y == 0;
cvx_end

mu = mu + rho*(x - y);
csvwrite(strcat(arg1,"y.csv"),y)
csvwrite(strcat(arg1,"mu.csv"),mu)
dlmwrite('observacion.dat',y','-append');

end
