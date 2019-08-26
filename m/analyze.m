clc
close all

path = '../data/25_08_19_5/';

R_line= [1.632,1.088,1.088];
X_line= [1.1019,.7346,.7346];

obs0 = csvread(strcat(path,'0/observacion.dat'));
obs1 = csvread(strcat(path,'1/observacion.dat'));
obs2 = csvread(strcat(path,'2/observacion.dat'));
obs3 = csvread(strcat(path,'3/observacion.dat'));

x0 = csvread(strcat(path,'0/estado.dat'));
x1 = csvread(strcat(path,'1/estado.dat'));
x2 = csvread(strcat(path,'2/estado.dat'));
x3 = csvread(strcat(path,'3/estado.dat'));

A0 = csvread(strcat(path,'0/A.csv'));
A1 = csvread(strcat(path,'1/A.csv'));
A2 = csvread(strcat(path,'2/A.csv'));
A3 = csvread(strcat(path,'3/A.csv'));

r0 = csvread(strcat(path,'0/residuo.dat'));
r1 = csvread(strcat(path,'1/residuo.dat'));
r2 = csvread(strcat(path,'2/residuo.dat'));
r3 = csvread(strcat(path,'3/residuo.dat'));

r_d0 = csvread(strcat(path,'0/residuo_dual.dat'));
r_d1 = csvread(strcat(path,'1/residuo_dual.dat'));
r_d2 = csvread(strcat(path,'2/residuo_dual.dat'));
r_d3 = csvread(strcat(path,'3/residuo_dual.dat'));

[v_ancestor0,v0,p0,q0,P0,Q0,l0,P10,Q10,l10] = load_node_root(strcat(path,'0/estado.dat'), 1, 200);
[v_ancestor1,v1,p1,q1,P1,Q1,l1,P11,Q11,l11,P21,Q21,l21] = load_node1(strcat(path,'1/estado.dat'), 1, 200);
[v_ancestor2,v2,p2,q2,P2,Q2,l2] = load_node_leaf(strcat(path,'2/estado.dat'), 1, 200);
[v_ancestor3,v3,p3,q3,P3,Q3,l3] = load_node_leaf(strcat(path,'3/estado.dat'), 1, 200);

x0 = [v_ancestor0,v0,p0,q0,P0,Q0,l0,P10,Q10,l10];
x1 = [v_ancestor1,v1,p1,q1,P1,Q1,l1,P11,Q11,l11,P21,Q21,l21];
x2 = [v_ancestor2,v2,p2,q2,P2,Q2,l2];
x3 = [v_ancestor3,v3,p3,q3,P3,Q3,l3];

[v_ancestor0,v0,p0,q0,P0,Q0,l0,P10,Q10,l10] = load_node_root(strcat(path,'0/observacion.dat'), 1, 200);
[v_ancestor1,v1,p1,q1,P1,Q1,l1,P11,Q11,l11,P21,Q21,l21] = load_node1(strcat(path,'1/observacion.dat'), 1, 200);
[v_ancestor2,v2,p2,q2,P2,Q2,l2] = load_node_leaf(strcat(path,'2/observacion.dat'), 1, 200);
[v_ancestor3,v3,p3,q3,P3,Q3,l3] = load_node_leaf(strcat(path,'3/observacion.dat'), 1, 200);


y0 = [v_ancestor0,v0,p0,q0,P0,Q0,l0,P10,Q10,l10];
y1 = [v_ancestor1,v1,p1,q1,P1,Q1,l1,P11,Q11,l11,P21,Q21,l21];
y2 = [v_ancestor2,v2,p2,q2,P2,Q2,l2];
y3 = [v_ancestor3,v3,p3,q3,P3,Q3,l3];

S_base = 100;

index = 100;
v_bus = [v0(index) v1(index) v2(index) v3(index)]
p_iny = S_base*[p0(index) p1(index) p2(index) p3(index)]
q_iny = S_base*[q0(index) q1(index) q2(index) q3(index)]

p_line = S_base* [P0(index) P1(index) P2(index) P3(index)]
q_line = S_base* [Q0(index) Q1(index) Q2(index) Q3(index)]
loss_line =  [0 R_line(1)*l1(index) R_line(2)*l2(index) R_line(3)*l3(index)]
current_line = [l0(index) l1(index) l2(index) l3(index)]

x0(index,:)
size(x0(index,:))
size(y0(index,:))
size(A0)
A3*x3(index,:)'
A3*y3(index,:)'

figure()
plot(x1(:,2:7))


% Residuo global
r= [r0,r1,r2,r3];
norma = zeros(index,1);
for i=1:1:index
    norma(i) = norm(r(i,:),2);
end


figure()
plot(norma)

P1(index)
p1(index)+P2(index)-l2(index)*R_line(2) +P3(index) - l3(index)*R_line(3)


