% #1 
clear;clc;
fprintf('Exercise #1\n\n');
omega = [0.5; 0.5; 0.7];
theta = 30/180*pi;

% Calculate rotation matrix using both methods:
fprintf('Calculating rotation matrix using expm:\n');
omegahat = ssm(omega);                  % See ssm.m
R = expm(omegahat*theta)    
fprintf('...and now using the Rodrigues formula:\n');
R = rodrigues(omega, theta)             % See rodrigues.m

fprintf('We can see that the expm function loses');  
fprintf(' some accuracy to rounding errors.\n');
fprintf('\n\nstrike a key to continue\n');
pause;

% #2
clear;clc;
fprintf('Exercise #2\n\n');
R = [0.6350 -0.0520 0.7708; 0.4676 0.8200 -0.3308; -0.6149 0.57 0.545]
fprintf('For rotation matrices, R*R'' = I\n')
R*R' 
fprintf('This is really close.\n')
fprintf('\nFor rotation matrices, det(R) = 1\n');
det(R)          
fprintf('This is also really close.\n');

fprintf('\n\nstrike a key to continue\n\n');
pause;
clc;

fprintf('Exercise #2 (cont.)\n\n');
fprintf('Now obtaining theta:\n');
theta = acos((trace(R) - 1)/2)          %theta = pi/3.  Cool.
fprintf('Solving for omegahat using logm:\n');
omegahat  = logm(R)/theta               
fprintf('...and now using our formula:\n');
omega = 1/(2*sin(theta))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)] 
fprintf('So they agree to within rounding errors.\n');

fprintf('\n\nstrike a key to continue\n');
pause;

% #3
clear;clc;
fprintf('Exercise #3\n\n');

fprintf('Execute the following rotations:\n')
fprintf('(a) Rotate around the x-axis by epsilon radians\n')
fprintf('(b) Rotate around the y-axis by epsilon radian\n')
fprintf('(c) Rotate around the x-axis by -epsilon radians\n')
fprintf('(d) Rotate around the y-axis by -epsilon radians\n')

for epsilon = [0.5 .1 .01 .001 .0001]  
    epsilon
    [R, w] = rotation(epsilon)         
    fprintf('\nstrike a key to continue\n\n');
    pause;                                    
    clc;
    fprintf('Exercise #3 (cont.)\n\n');
end;                                     
fprintf('As epsilon goes to zero, the rotation matrix R'); 
fprintf(' approaches I (no rotation),\nand the axis of'); 
fprintf(' rotation w approaches a unit vector in the negative z\n'); 
fprintf('direction [0;0;-1] (see rotation.m)\n');

fprintf('\n\nstrike a key to continue\n');
pause;                                  
                                       
% #4 (See screw.m)

% #5 
clear;clc;
fprintf('Exercise #5\n\n');
w = [0.5; 0.5; 0.7]
v = [0.3; 0.8; 0.5196]
fprintf('Calculating the corresponding homogeneous transformation:\n')
g = screw(w,v)