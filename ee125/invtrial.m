clear; clc;

w1=[1/sqrt(8);sqrt(3/8);1/sqrt(2)];
w2=[-1/sqrt(3);0;sqrt(2/3)];
r=[1;2;-3];
% (w1 and w2 both pass through r).

p=[5;-7;12];
q1=[17.3037;-3.1128;2.48175];
q2=[14.6711;-9.34914;-0.490328];


[theta1, theta2] = sp2(w1,w2,p-r,q2-r)
theta = sp1(w1,p-r,q1-r)

w1 = [0;1;0];
w2 = [0;1;0];
p = [10;0;0];
q = [-5;0;5];
r1 = [0;100;0];
r2 = [5;-100;0];
[theta1 theta2] = sp2prime(w1,w2,p,q,r1,r2)