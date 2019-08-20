function [] = update_y(arg1, arg2)

%fprintf(strcat(arg1,"A.csv"))

A = csvread(strcat(arg1,"A.csv"));
x = csvread(strcat(arg1,"x.csv"));
y = csvread(strcat(arg1,"y.csv"));
mu = csvread(strcat(arg1,"mu.csv"));

B=[1,0,0]';

%%% Init variables OPF
n = size(y);
rho=arg2;

% Problem OPF: 
cvx_begin quiet
    variable y(n,1)
    %minimize(-mu'*y  + .5*rho*sum_square(x-y));
    minimize(-mu'*y  + .5*rho*(x-y)'*(x-y));
    A*y + B == 0;
cvx_end

mu = mu + rho*(x - y);
csvwrite(strcat(arg1,"y.csv"),y)
csvwrite(strcat(arg1,"mu.csv"),mu)
dlmwrite('observacion.dat',y','-append');
dlmwrite('multiplicadores.dat',mu','-append');

end
