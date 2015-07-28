function corio = corio(theta, thetadot, param)

l1 = param(1);
l2 = param(2);
m1 = param(3);
m2 = param(4);
t1d = thetadot(1);
t2d = thetadot(2);

corio =  [-(l1*l2*m2*t2d*sin(theta(2)))/2,...
    -(l1*l2*m2*(t1d + t2d)*sin(theta(2)))/2;... 
    (l1*l2*m2*t1d*sin(theta(2)))/2, 0];
