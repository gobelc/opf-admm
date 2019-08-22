function [] = update_x(arg1,arg2)

%fprintf(strcat(arg1,"A.csv"))
A = csvread(strcat(arg1,"A.csv"));
A = A(1:2,:);
x_old = csvread(strcat(arg1,"x.csv"));
y = csvread(strcat(arg1,"y.csv"));
mu = csvread(strcat(arg1,"mu.csv"));

%%% Init variables OPF
n = size(y,1);
c = zeros(n,1);
c(3)=1;

rho=arg2

n_childs = floor((n-7)/3);

% Problem OPF: 
cvx_begin
    variable x(n,1)
    minimize(square(c'*x) + mu'*x + .5*rho*(x-y)'*(x-y));
    %x(1)==x(2);
    x(1)==1.
    x(2)==1.;
    %x(3)<=10;
    %x(3)>=0.;
    %x(5)==0.;
    %x(6)==0.;
    %x(7)==0.;
cvx_end
csvwrite(strcat(arg1,"x.csv"),x)
residuo = norm(x-y,2);
dlmwrite('estado.dat',x','-append');
dlmwrite('residuo.dat',residuo,'-append');
dlmwrite('costo.dat',c'*x,'-append');
end
