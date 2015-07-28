clear;clc;

syms L0 L1 real;
syms th real;

omega = [-1 0 0]';
q2 = [0 0 0]';
q3 = [0 L0 0]';
p2 = [q2; 1];

g_st0 = [1 0 0 0; 0 1 0 L0+L1; 0 0 1 0; 0 0 0 1];
v = cross(-omega, q3);
g_3 = screw(omega, v, -th);

simple(inv(g_st0)*g_3*p2)