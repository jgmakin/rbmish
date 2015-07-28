function potfor = potfor(theta,param)

l1 = param(1);
l2 = param(2);
m1 = param(3);
m2 = param(4);
grav = param(5);
t1 = theta(1);
t2 = theta(2);
c1 = cos(t1);
c12 = cos(t1 + t2);

potfor =   [grav*((l1*m1*c1)/2 + l1*m2*c1 + (l2*m2*c12)/2), 
   (grav*l2*m2*c12)/2];
