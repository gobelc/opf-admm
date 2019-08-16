clc
close all

[v_ancestor0,v0,p0,q0,P10,Q10,l10,P20,Q20,l20] = load_node_root('../data/16_08_19_2/0/estado.dat', 1, 200);
[v_ancestor1,v1,p1,q1,P11,Q11,l11,P21,Q21,l21] = load_node1('../data/16_08_19_2/1/estado.dat', 1, 200);
[v_ancestor2,v2,p2,q2,P12,Q12,l12] = load_node_leaf('../data/16_08_19_2/2/estado.dat', 1, 200);
[v_ancestor3,v3,p3,q3,P13,Q13,l13] = load_node_leaf('../data/16_08_19_2/3/estado.dat', 1, 200);


v_bus = [v0(200) v1(200) v2(200) v3(200)]
p_iny = [p0(200) p1(200) p2(200) p3(200)]
q_iny = [q0(200) q1(200) q2(200) q3(200)]