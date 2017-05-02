function R = assembleLissajousTensor(x0,v0,xmin,xmax,T,dt,eta)
% R = assembleLissajousTensor(x0,v0,smin,smax,Ncases,T,dt,eta)
%
% Given matrices (Ncases x Ndims) for the initial positions (x0) and
% velocities (v0); the boundaries of the space (xmin,xmax); the number of
% time steps (T), and the discrete interval (dt); generate a series of
% Ncases Lissajous curves, each of length T.  

% The fudge-factor eta is for margin control: it is the percentage of the 
% space that is actually occupied by the curves.
%
% NB: the reference trajectories live in 90% of joint-angle space; this 
% gives the noisy outputs a little breathing room wrt the margins.

%-------------------------------------------------------------------------%
% Revised: 01/08/14
%   -brought into its own function
%   -eliminated Ncases as in input argument
% Created: 12/09/13
%   by JGM
%-------------------------------------------------------------------------%


%%%%%%
%%% arcane: you should really check to make sure the ICs are never drawn
%%% from points outside 90% (eta) of the response area
%%%%%%


% params
Ncases = size(x0,1);
Ninputs = size(xmin,1);


% Lissajous parameters
t = 0:dt:(dt*(T-1));                            % 1 x T
b = (xmax' + xmin')/2;                          % 1 x Ninputs
A = (xmax' - xmin')/2;                          % 1 x Ninputs
A = eta*A;                                      % ...margin ctrl
z0 = (x0-b)./repmat(A,[Ncases,1]);              % Ncases x Ninputs
phi = asin(z0);                                 % Ncases x Ninputs
w = (asin(dt*(v0./A) + z0)-phi)/dt;             % Ncases x Ninputs

R = A.*sin(reshape(w(:)*t,[Ncases,Ninputs,T]) + phi) + b;


end





% %-------------------------------------------------------------------------%
% function R = assembleReferenceTensor2(wmax,smin,smax,Ncases,T,dt)
% % This version uses random w and phi rather than random x0 and v0
% 
% % fixed param
% Ninputs = size(smin,1);
% dmax = 2*pi;
% eta = 1.2;      % fudge factor for margin control
% 
% % get (possibly random) parameters
% phi = dmax*rand(Ncases,Ninputs);          % Ncases x 1 (see [1])
% w = wmax*rand(Ncases,Ninputs); % 7/3, 10/3  % Ncases x Ninputs
% t = 0:dt:(dt*(T-1));                        % 1 x T
% A = ones(1,Ninputs);                        % 1 x Ninputs (ones -> rand?)
% 
% % put into tensors that you can multiply properly
% pTensor = repmat(phi,[1,1,T]);
% ATensor = repmat(A,[Ncases,1,T]);
% wTensor = repmat(w,[1,1,T]);
% tTensor = repmat(shiftdim(t,-1),[Ncases,Ninputs]);
% 
% unitR = ATensor.*sin(wTensor.*tTensor + pTensor);
% 
% R = scalefxn(longdata(unitR)',-eta*ones(Ninputs,1),eta*ones(Ninputs,1),smin,smax)';
% R = shortdata(Ncases,3,R);
% 
% end
% %-------------------------------------------------------------------------%