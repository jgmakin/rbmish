function dynamics = setDynamics(Ndims,mods,xmin,xmax)
% setDynamics   Set dynamics for rEFH simulations
%
% NB that this function sets a *uniform* prior over position and a tight
% prior about 0 over velocity.

%-------------------------------------------------------------------------%
% Revised: 08/24/16
%   -changed input arguments from params to the necessary parts
%   -commented
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%


% init
NnonECpops = sum(~strcmp(mods,'Efference-Copy'));
Nstates = 2*Ndims;                              % second-order dynamics
xrange = xmax - xmin;


% x[t+1] = A*x[t] + B*noise + muX,      y[t] = C*x[t] + noise
dt = 0.05;
m = 5;
b = 0.25;                                       % damping
k = 3;                                          % spring force 
A = [eye(Ndims)             dt*eye(Ndims);...
    -k/m*dt*eye(Ndims)      -(b/m*dt-1)*eye(Ndims)];
C = zeros(Ndims*NnonECpops,Nstates);
C(1:Ndims*NnonECpops,1:Ndims*NnonECpops) = eye(Ndims*NnonECpops);

% transition noise
posVar  = 5e-7;
velVar  = 5e-5;
SigmaX  = blkdiag(posVar*eye(Ndims),velVar*eye(Ndims));

% prior---assumes no position-velocity coupling
muX0    = xrange/2 + xmin;
SigmaX0 = diag(Inf*(1:Ndims)); % diag(range/4);
muV0    = zeros(Ndims,1);
SigmaV0 = 5e-10;


% are there controls?
if any(strcmp(mods,'Efference-Copy'));
    
    % x[t+1] = A*x[t] + G*u[t] + noise,     y[t] = C*x[t] + noise 
    % u[t+1] = F*u[t] + noise,              v[t+1] = H*u[t] + noise
    
    Ninputs = Ndims;
    alpha = 0.9994;
    F = alpha*eye(Ninputs);
    G = zeros(Nstates,Ninputs);
    G((end-Ninputs+1):end,(end-Ninputs+1):end) = eye(Ninputs)*dt/m;
    H = eye(Ninputs);
    
    muU0 = 0*ones(Ninputs,1);
    SigmaU0 = (5e-10)*eye(Ninputs);

    % put position/velocity noise in range that yields good HSVs
    SigmaU = 7.5e-4*eye(Ninputs);
    SigmaX(1:Ndims,1:Ndims) = 5e-5*eye(Ndims);
    SigmaX(Ndims+1:end,Ndims+1:end) = 1e-6*eye(Ndims);
    
    % load into structure
    dynamics.muU0 = muU0;
    dynamics.SigmaU0 = SigmaU0;
    
    % augment the parameters
    A = [A, G; zeros(size(F,1),size(A,2)), F];
    C = blkdiag(C,H);
    SigmaX = blkdiag(SigmaX,SigmaU);
    
end

% test Hankel singular values (can the system be model-order reduced?)
try
    P = dlyap(A,SigmaX);                         % SigmaX = B*B'
    Q = dlyap(A',C'*C);
    hsvs = sqrt(eig(P*Q));
    fprintf('Hankel singular values: %f and %f\n',hsvs)
catch ME
    fprintf('can''t print Hankel SVs: out of control ');
    fprintf('toolbox licenses -- jgm\n')
end

% finally, create state-transition "bias" at all zeros
muX = zeros(size(SigmaX,1),1);

% load into a structure
dynamics.A = A;
dynamics.C = C;
dynamics.SigmaX = SigmaX;
dynamics.muX = muX;
dynamics.muX0 = muX0;
dynamics.SigmaX0 = SigmaX0;
dynamics.SigmaV0 = SigmaV0;
dynamics.muV0 = muV0;
dynamics.dt = dt;
dynamics.m = m;




end


% *** (1) ***
% A RECIPE FOR SETTING "INTERESTING" DYNAMICS
%
% (0) pick a sampling interval, dt, for which the DT approx. is good
%   (*) I chose 0.05 s.
%
% (1) pick the efference-copy gain to make it compatible with dt.
%   (*) I chose 10: any lower and you'll start to get examples w/no spikes.
%   (a) This yields 10 spikes/0.05 s = 200 spikes/s, which is reasonable.
%   (b) It gives covariances of about ...
%
% (2) pick muU0 to be zero
%   (a) this will prevent the system from running away....
% (3) pick SigmaU0 to be small (compared to the efference-copy precision)
%   (a) you don't want initial position to matter much: you want to test
%   how well the system(s) learn dynamics, not the p0(u,x).
%   (*) you chose 1e-10 (cf. 1b).
%
% (4) set control "decay," alpha, so that the dynamics are significant
%   (a) e.g., you have to prevent u from decaying to zero
%   (*) you chose 1, i.e. no decay
%
% (5) set the control transition noise, SigmaU, so that:
%   (a) the total motion of u is w/i, or close to, the decoding precision
%       => useful to track the *dynamics* of u
%       => small SigmaU
%   (b) u doesn't just stay "close" to zero
%       => useful to assume a control exists
%       => large SigmaU
%   (c) HOW CLOSE IS "close"?
%   (*) you chose 5e-8

%   (b) u doesn't just stay close to its IC
%       => SigmaU sufficiently high
%   (c) the total movement of u is within the variance of uhat; this will
%   ensure that it is useful 
























