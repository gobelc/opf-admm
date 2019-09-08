function [] = update_x(arg1,arg2)

%fprintf(strcat(arg1,"A.csv"))
A = csvread(strcat(arg1,"A.csv"));
A = A(1:2,:);
x_old = csvread(strcat(arg1,"x.csv"));
y = csvread(strcat(arg1,"y.csv"));
mu = csvread(strcat(arg1,"mu.csv"));

%%% Init variables OPF
S_base = 100; %MVA
n = size(y,1);
c1 = zeros(n,1);
c2 = zeros(n,1);

%% COST FUNCTION 
R = 0.00304;
%c(7) = R; %Minimize power loss
c1(3) = 1; %Minimize power loss
c2(4) = 1; %Minimize power loss

%Restricciones de linea
Pmax=100;
Qmax=100;

rho=arg2

n_childs = floor((n-7)/3);

% Problem OPF: 
cvx_begin
    variable x(n,1)
    minimize(norm(c1'*x)+norm(c2'*x) + mu'*x + .5*rho*(x-y)'*(x-y));
    x(1)==1.;
    x(2)==1.;
    x(3)<=100/S_base;
    x(5)==0.;
    x(6)==0.;
    x(7)==0.;
    x(8)==y(8);
    x(9)==y(9);
    x(10)==y(10);    
cvx_end
costo = norm(c1'*x)+norm(c2'*x);
csvwrite(strcat(arg1,"x.csv"),x)
residuo = norm(x-y,2);
dlmwrite('estado.dat',x','-append');
dlmwrite('residuo.dat',residuo,'-append');
dlmwrite('costo.dat',costo,'-append');
end
