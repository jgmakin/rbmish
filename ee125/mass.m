function mass = mass(theta, param)

l1 = param(1);
l2 = param(2);
m1 = param(3);
m2 = param(4);

mass = [(l2^2*m2)/4 + l1^2*(m1/4 + m2) + l1*l2*m2*cos(theta(2)), ...
    (l2^2*m2)/4 + (l1*l2*m2*cos(theta(2)))/2;.... 
    (l2^2*m2)/4 + (l1*l2*m2*cos(theta(2)))/2, (l2^2*m2)/4];