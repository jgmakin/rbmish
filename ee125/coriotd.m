function coriotd = coriotd(theta, thetadot, param)

l1 = param(1);
l2 = param(2);
m1 = param(3);
m2 = param(4);
t1d = thetadot(1);
t2d = thetadot(2);

coriotd = [-(l1*l2*m2*t2d*(2*t1d + t2d)*sin(theta(2)))/2;...
          (l1*l2*m2*t1d^2*sin(theta(2)))/2];
 