clc
close all

path = '../data/19_08_19_4/';

A0 = csvread(strcat(path,'0/A.csv'));
A1 = csvread(strcat(path,'1/A.csv'));
A2 = csvread(strcat(path,'2/A.csv'));
A3 = csvread(strcat(path,'3/A.csv'));
r0 = csvread(strcat(path,'0/residuo.dat'));
r1 = csvread(strcat(path,'1/residuo.dat'));
r2 = csvread(strcat(path,'2/residuo.dat'));
r3 = csvread(strcat(path,'3/residuo.dat'));
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


index = 120;
v_bus = [v0(index) v1(index) v2(index) v3(index)]
p_iny = [p0(index) p1(index) p2(index) p3(index)]
q_iny = [q0(index) q1(index) q2(index) q3(index)]

x0(index,:)
size(x0(index,:))
size(y0(index,:))
size(A0)
A3*x3(index,:)'
A3*y3(index,:)'

figure()
plot(x1(:,2:7))
