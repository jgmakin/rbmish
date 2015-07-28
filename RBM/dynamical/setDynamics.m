function dynamics = setDynamics(params)

% init
NnonECpops = sum(~strcmp(params.mods,'Efference-Copy'));
Ndims = params.Ndims;
Nstates = 2*Ndims;                              % second-order dynamics
xmin = params.smin(:,strcmp(params.mods,params.NS));
xmax = params.smax(:,strcmp(params.mods,params.NS));
xrange = xmax - xmin;


% x[t+1] = A*x[t] + B*noise + muX,      y[t] = C*x[t] + noise
dt = 0.05;
m = 5;
b = 0.25;                                       % damping
k = 3;                                          % spring force 
A = [eye(Ndims)             dt*eye(Ndims);...
    -k/m*dt*eye(Ndims)      -(b/m*dt-1)*eye(Ndims)];
%%%%%%
% A = [0.9892 0.0474; -0.3361 0.8269];
% A = [1 0.05; -0.34 0.83];
%%%%%%
C = zeros(Ndims*NnonECpops,size(A,1));
C(1:Ndims*NnonECpops,1:Ndims*NnonECpops) = eye(Ndims*NnonECpops);

% transition noise
posVar = 5e-7;
velVar = 5e-5;
SigmaX = [posVar*eye(Ndims) zeros(Ndims); zeros(Ndims) velVar*eye(Ndims)];
muX = zeros(Nstates,1);

% prior---assumes no position-velocity coupling
muX0 = xrange/2 + xmin;
SigmaX0 = diag(Inf*(1:Ndims)); % diag(range/4);
muV0 = zeros(Ndims,1);
SigmaV0 = 5e-10;


% are there controls?
if any(strcmp(params.mods,'Efference-Copy'));
    
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

    % put noise in range that yields good HSVs
    SigmaU = 7.5e-4*eye(Ninputs);
    SigmaX(1:Ndims,1:Ndims) = 5e-5*eye(Ndims);
    SigmaX(Ndims+1:end,Ndims+1:end) = 1e-6*eye(Ndims);
    
    
    
    % load into structure
    dynamics.F = F;
    dynamics.G = G;
    dynamics.H = H;
    dynamics.SigmaU = SigmaU;
    dynamics.muU0 = muU0;
    dynamics.SigmaU0 = SigmaU0;

    
    % test Hankel singular values (can the system be model-order reduced?)
    Gamma = zeros(Nstates+Ninputs,Nstates+Ninputs);
    Gamma(1:Nstates,1:Nstates) = A; 
    Gamma(1:Nstates,Nstates+1:end) = G;
    Gamma(Nstates+1:end,Nstates+1:end) = F;
    %%% Lambda = zeros(2*Ndims,Nstates+Ninputs);
    %%% Lambda(1:Ndims,1:Nstates) = C;
    %%% Lambda(Ndims+1:end,Nstates+1:end) = H; % we don't care how visible
    %%%the control is in *its* observation, which we're sure is big
    Lambda = zeros(Ndims*NnonECpops,Nstates+Ninputs);
    Lambda(1:Ndims*NnonECpops,1:Nstates) = C;  
    %%%
    UpsilonZ = zeros(Nstates+Ninputs,Nstates+Ninputs);
    UpsilonZ(1:Nstates,1:Nstates) = SigmaX;
    UpsilonZ(Nstates+1:end,Nstates+1:end) = SigmaU;
    try 
        P = dlyap(Gamma,UpsilonZ);                  % UpsilonZ = B*B'
        Q = dlyap(Gamma',Lambda'*Lambda);
        hsvs = sqrt(eig(P*Q));
        fprintf('Hankel singular values: %f, %f, and %f\n',hsvs)
    catch ME
        fprintf('can''t print Hankel SVs: out of control ');
        fprintf('toolbox licenses -- jgm\n')
    end
    
else
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
end


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

% what happens at the walls?
% bouncing | resetting | linearInXY | quadraticbowl | repelling | wrapping
dynamics.walls = 'wrapping';





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
























