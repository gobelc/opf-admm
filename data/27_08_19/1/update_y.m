function [] = update_y(arg1,arg2)

%fprintf(strcat(arg1,"A.csv"))

A = csvread(strcat(arg1,"A.csv"));
x = csvread(strcat(arg1,"x.csv"));
y = csvread(strcat(arg1,"y.csv"));
mu = csvread(strcat(arg1,"mu.csv"));


%%% Init variables OPF
n = size(y);
rho=arg2;
y_old = y;


% Problem OPF: 
cvx_begin quiet
    variable y(n,1)
    minimize(-mu'*y  + .5*rho*(x-y)'*(x-y));
    A*y == 0;
cvx_end

mu = mu + rho*(x - y);
residuo = norm(x-y);
residuo_dual = rho * norm(y-y_old,2);
tau = 2;

if residuo> 10* residuo_dual
    rho = tau*rho;
elseif residuo_dual> 10 * residuo
    rho = rho / tau;
end       
    
csvwrite(strcat(arg1,"y.csv"),y)
csvwrite(strcat(arg1,"mu.csv"),mu)
dlmwrite('rho.csv',rho);
dlmwrite('residuo_dual.dat',residuo_dual,'-append');
dlmwrite('observacion.dat',y','-append');
dlmwrite('multiplicadores.dat',mu','-append');
end
