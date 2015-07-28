clear;clc;
syms L0 L1 th1 th2 th3;


w1 = [0 0 1]';
v1 = [0 0 0]';
w2 = [-1 0 0]';
v2 = [0 -L0 0]';
w3 = [0 0 0]';
v3 = [0 1 0]';

one = screw(w1,v1,th1);
two = screw(w2,v2,th2);
three = screw(w3,v3,th3);
four = [1 0 0 0; 0 1 0 L1; 0 0 1 L0; 0 0 0 1]
one*two*three*four

% M1 = dh(0,L0,pi/2,th1);
% M2 = dh(0,0,pi/2,th2);
% M3 = dh(0,L1,pi/2,th3);
% M = M1*M2*M3