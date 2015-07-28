function xdot = robotsim(t,x)

roboparam = [4.0;3.0;2.0;3.5;9.81];

minv = inv(mass(x, roboparam));
ctd = coriotd(x, [x(3);x(4)], roboparam);
pf = potfor(x, roboparam);
torque = [x(5);x(6)];

xdot = [x(3);x(4);minv * (torque- ctd - pf);0;0];
