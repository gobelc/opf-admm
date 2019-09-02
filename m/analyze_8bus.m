clc
close all
index = 100;

path = '../data/1_09_19_2/';

R_line= [0.00304,0.00304,0.00304,0.00304,0.00304,0.00304,0.00304];
X_line= [0.0304,0.0304,0.0304,0.0304,0.0304,0.0304,0.0304];

obs0 = csvread(strcat(path,'0/observacion.dat'));
obs1 = csvread(strcat(path,'1/observacion.dat'));
obs2 = csvread(strcat(path,'2/observacion.dat'));
obs3 = csvread(strcat(path,'3/observacion.dat'));
obs4 = csvread(strcat(path,'4/observacion.dat'));
obs5 = csvread(strcat(path,'5/observacion.dat'));
obs6 = csvread(strcat(path,'6/observacion.dat'));
obs7 = csvread(strcat(path,'7/observacion.dat'));

x0 = csvread(strcat(path,'0/estado.dat'));
x1 = csvread(strcat(path,'1/estado.dat'));
x2 = csvread(strcat(path,'2/estado.dat'));
x3 = csvread(strcat(path,'3/estado.dat'));
x4 = csvread(strcat(path,'4/estado.dat'));
x5 = csvread(strcat(path,'5/estado.dat'));
x6 = csvread(strcat(path,'6/estado.dat'));
x7 = csvread(strcat(path,'7/estado.dat'));

A0 = csvread(strcat(path,'0/A.csv'));
A1 = csvread(strcat(path,'1/A.csv'));
A2 = csvread(strcat(path,'2/A.csv'));
A3 = csvread(strcat(path,'3/A.csv'));
A4 = csvread(strcat(path,'4/A.csv'));
A5 = csvread(strcat(path,'5/A.csv'));
A6 = csvread(strcat(path,'6/A.csv'));
A7 = csvread(strcat(path,'7/A.csv'));

r0 = csvread(strcat(path,'0/residuo.dat'));
r1 = csvread(strcat(path,'1/residuo.dat'));
r2 = csvread(strcat(path,'2/residuo.dat'));
r3 = csvread(strcat(path,'3/residuo.dat'));
r4 = csvread(strcat(path,'4/residuo.dat'));
r5 = csvread(strcat(path,'5/residuo.dat'));
r6 = csvread(strcat(path,'6/residuo.dat'));
r7 = csvread(strcat(path,'7/residuo.dat'));

r_d0 = csvread(strcat(path,'0/residuo_dual.dat'));
r_d1 = csvread(strcat(path,'1/residuo_dual.dat'));
r_d2 = csvread(strcat(path,'2/residuo_dual.dat'));
r_d3 = csvread(strcat(path,'3/residuo_dual.dat'));
r_d4 = csvread(strcat(path,'4/residuo_dual.dat'));
r_d5 = csvread(strcat(path,'5/residuo_dual.dat'));
r_d6 = csvread(strcat(path,'6/residuo_dual.dat'));
r_d7 = csvread(strcat(path,'7/residuo_dual.dat'));

[v_ancestor0,v0,p0,q0,P0,Q0,l0,P10,Q10,l10] = load_node_root(strcat(path,'0/estado.dat'), 1, index);
[v_ancestor1,v1,p1,q1,P1,Q1,l1,P11,Q11,l11,P21,Q21,l21] = load_node_hub(strcat(path,'1/estado.dat'), 1, index);
[v_ancestor2,v2,p2,q2,P2,Q2,l2,P12,Q12,l12] = load_node_bridge(strcat(path,'2/estado.dat'), 1, index);
[v_ancestor3,v3,p3,q3,P3,Q3,l3,P13,Q13,l13,P23,Q23,l23] = load_node_hub(strcat(path,'3/estado.dat'), 1, index);
[v_ancestor4,v4,p4,q4,P4,Q4,l4,P14,Q14,l14] = load_node_bridge(strcat(path,'4/estado.dat'), 1, index);
[v_ancestor5,v5,p5,q5,P5,Q5,l5] = load_node_leaf(strcat(path,'5/estado.dat'), 1, index);
[v_ancestor6,v6,p6,q6,P6,Q6,l6] = load_node_leaf(strcat(path,'6/estado.dat'), 1, index);
[v_ancestor7,v7,p7,q7,P7,Q7,l7] = load_node_leaf(strcat(path,'7/estado.dat'), 1, index);

x0 = [v_ancestor0,v0,p0,q0,P0,Q0,l0,P10,Q10,l10];
x1 = [v_ancestor1,v1,p1,q1,P1,Q1,l1,P11,Q11,l11,P21,Q21,l21];
x2 = [v_ancestor2,v2,p2,q2,P2,Q2,l2,P12,Q12,l12];
x3 = [v_ancestor3,v3,p3,q3,P3,Q3,l3,P13,Q13,l13,P23,Q23,l23];
x4 = [v_ancestor4,v4,p4,q4,P4,Q4,l4,P14,Q14,l14];
x5 = [v_ancestor5,v5,p5,q5,P5,Q5,l5];
x6 = [v_ancestor6,v6,p6,q6,P6,Q6,l6];
x7 = [v_ancestor7,v7,p7,q7,P7,Q7,l7];

[v_ancestor0,v0,p0,q0,P0,Q0,l0,P10,Q10,l10] = load_node_root(strcat(path,'0/observacion.dat'), 1, index);
[v_ancestor1,v1,p1,q1,P1,Q1,l1,P11,Q11,l11,P21,Q21,l21] = load_node_hub(strcat(path,'1/observacion.dat'), 1, index);
[v_ancestor2,v2,p2,q2,P2,Q2,l2,P12,Q12,l12] = load_node_bridge(strcat(path,'2/observacion.dat'), 1, index);
[v_ancestor3,v3,p3,q3,P3,Q3,l3,P13,Q13,l13,P23,Q23,l23] = load_node_hub(strcat(path,'3/observacion.dat'), 1, index);
[v_ancestor4,v4,p4,q4,P4,Q4,l4,P14,Q14,l14] = load_node_bridge(strcat(path,'4/observacion.dat'), 1, index);
[v_ancestor5,v5,p5,q5,P5,Q5,l5] = load_node_leaf(strcat(path,'5/observacion.dat'), 1, index);
[v_ancestor6,v6,p6,q6,P6,Q6,l6] = load_node_leaf(strcat(path,'6/observacion.dat'), 1, index);
[v_ancestor7,v7,p7,q7,P7,Q7,l7] = load_node_leaf(strcat(path,'7/observacion.dat'), 1, index);

y0 = [v_ancestor0,v0,p0,q0,P0,Q0,l0,P10,Q10,l10];
y1 = [v_ancestor1,v1,p1,q1,P1,Q1,l1,P11,Q11,l11,P21,Q21,l21];
y2 = [v_ancestor2,v2,p2,q2,P2,Q2,l2,P12,Q12,l12];
y3 = [v_ancestor3,v3,p3,q3,P3,Q3,l3,P13,Q13,l13,P23,Q23,l23];
y4 = [v_ancestor4,v4,p4,q4,P4,Q4,l4,P14,Q14,l14];
y5 = [v_ancestor5,v5,p5,q5,P5,Q5,l5];
y6 = [v_ancestor6,v6,p6,q6,P6,Q6,l6];
y7 = [v_ancestor7,v7,p7,q7,P7,Q7,l7];

S_base = 100;

v_bus = [v0(index) v1(index) v2(index) v3(index) v4(index) v5(index) v6(index) v7(index)]
p_iny = S_base*[p0(index) p1(index) p2(index) p3(index) p4(index) p5(index) p6(index) p7(index)]
q_iny = S_base*[q0(index) q1(index) q2(index) q3(index) q4(index) q5(index) q6(index) q7(index)]

p_line = S_base* [P0(index) P1(index) P2(index) P3(index) P4(index) P5(index) P6(index) P7(index)]
q_line = S_base* [Q0(index) Q1(index) Q2(index) Q3(index) Q4(index) Q5(index) Q6(index) Q7(index)]
loss_line =  [0 R_line(1)*l1(index) R_line(2)*l2(index) R_line(3)*l3(index) R_line(4)*l4(index) R_line(5)*l5(index) R_line(6)*l6(index) R_line(7)*l7(index)]
current_line = [l0(index) l1(index) l2(index) l3(index) l4(index) l5(index) l6(index) l7(index)]

figure()
plot(x3(:,2:7))


% Residuo global
r_primal= [r0(1:index),r1(1:index),r2(1:index),r3(1:index),r4(1:index),r5(1:index),r6(1:index),r7(1:index)];
r_dual = [r_d0(1:index),r_d1(1:index),r_d2(1:index),r_d3(1:index),r_d4(1:index),r_d5(1:index),r_d6(1:index),r_d7(1:index)];
norma_primal = zeros(index,1);
norma_dual = zeros(index,1);

for i=1:1:index
    norma_primal(i) = norm(r_primal(i,:),2);
    norma_dual(i) = norm(r_dual(i,:),2);
end


figure()
plot(norma_primal)
hold on;
plot(norma_dual)
grid on
legend('Primal','Dual')
title('Norma de residuo del primal')


%%% Chequear consistencia

%nodo 3
S_base=100;
%Potencia activa
A_p=[-2.95237e-38,0,1,0,-1,0,0,1,0,-0.00304,1,0,-0.00304];
fprintf('Residuo potencia activa: %f\n',S_base*A_p*y3(index,:)');
fprintf('Potencia inyectada: %f\n',S_base*y3(index,3));
fprintf('Potencia activa linea: %f\n',S_base*y3(index,5));
fprintf('Potencia activa hijos: %f\n',S_base*y3(index,8)+S_base*y3(index,11));
fprintf('Perdidas linea: %f\n',S_base*(A_p(10)*y3(index,10) + A_p(13)*y3(index,13)));