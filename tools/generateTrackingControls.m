function [Z,U] = generateTrackingControls(R,umin,umax,dynamics)
% [Z,U] = generateTrackingControls(R,umin,umax,dynamics)
%
%   This function finds a set of control, U, that force the system in the
%   structure dynamics to track the reference R, using a combination of
%   exact tracking (assuming the LTI system has relative degree two) and 
%   LQR control.  It then simulates that control system, using the noise 
%   properties also contained in dynamics, and returns the resulting
%   states, Z, in a tensor of size Ncases x Nstates x T.  Note that R must
%   have dimensions Ncases x Noutputs x T.


%-------------------------------------------------------------------------%
% Revised: 01/08/13
%   -replaced LQI (?) with LQR for the error (y-r) + exact control of the
%   Lissajous curves.  (The latter appears superfluous---shades of your
%   dissertation....)
%   -implemented exact solution of Lissajous tracking for any linear system
%   with relative degree 2 (I think...)
%   -broke out the LissajousCurve tensor generator into another file
%   -reduced input arguments so that this function is more general than
%   your RBM code
% Revised: 01/07/13
%   -added rails to the control.  This suggests perhaps using
%   model-predictive control...
% Revised: 01/06/13
%   -now gets "B" from params.dynamics.G rather than computing it here
% Revised: 12/18/13
%   -massively revised:
%       --removed hard-coding
%       --extended to all Ncases
%       --added noise
% Created: 12/09/13
%   by JGM
%-------------------------------------------------------------------------%


% the dynamical control system
A = dynamics.A;
B = dynamics.G;
C = dynamics.C;
SigmaX = dynamics.SigmaX;
Nstates = size(A,1);
Ninputs = size(B,2);
[Ncases,Noutputs,T] = size(R);
% D = zeros(Noutputs,Ninputs);


% pre-assemble state noise [this does indeed save time]
longW = mvnrnd(zeros(Ncases*T,size(SigmaX,1)),SigmaX);
W = shortdata(Ncases,3,longW); 
clear longW;


% control: get exact control/ICs this reference (ASSUMES RELDEG=2)
[Uopt,Zstar] = ExactTrackingControls(A,B,C,R);
%%% you may want to know this:
% 3*sqrt(var(longdata(Uopt)))

% control: linear-quadratic regulator for error induced by noise
QQ = 1*eye(Noutputs); %%% I just made up 6, I think
RR = 1*eye(Ninputs);
% [K,S,E] = dlqry(A,B,C,D,QQ,RR);
% [K,S,E] = dlqr(A,B,C'*QQ*C,RR);
if isa(A,'gpuArray')
    A = gather(A);
    P = dareJGM(A,B,RR,C'*QQ*C);
    P = gpuArray(P);
else
    P = dareJGM(A,B,RR,C'*QQ*C);
end

K = inv(RR)*B'*P; % K = RR\(B'*P);


% malloc
Z = zeros(Ncases,Nstates,T,'like',R);
U = zeros(Ncases,Ninputs,T,'like',R);

% init with the ideal IC
z = Zstar(:,:,1)';

% cycle through time
for t = 1:T
    
    % grab noise; compute error, input, and output
    w = W(:,:,t)';
    e = z - Zstar(:,:,t)';
    u = Uopt(:,:,t)' - K*e;
    
    % don't let the input go outside its limits
    u = (u >= umin).*u + (u < umin).*umin;
    u = (u <= umax).*u + (u > umax).*umax;
    
    % store input and state
    U(:,:,t) = u';
    Z(:,:,t) = z';

    % update state
    z = A*z + B*u + w;
 
end


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [Ustar,Zstar] = ExactTrackingControls(A,B,C,R)
% NB: THIS FUNCTION ASSUMES THAT THE RELATIVE DEGREE OF THE SYSTEM IS 2.
% It is, but you may not know that in general.
% 
% If you "knew" that the C matrix just pulls out the first two states, i.e.
% that (z1,z2) = (x,y); and that the A matrix forces (z3,z4) = (v_x,v_y);
% then you could compute the controls and ICs very simply in one step:
%
%   U = cat(3,[diff(R,2,3)/dt],zeros(Ncases,Ninputs,2));
%   z = [x0 v0]';
%
% But you will assume that you don't know this (it may change in the
% future, after all.  ((Notice that it would kind of make more sense to
% divide U by dt again, since the diff is taken twice: then it would be in
% units of acceleration.  You could then have multiplied B by dt to cancel
% this term---but anyway you didn't.))

% Assuming relative degree 2, the following constraints allow you to solve
% the control problem:
%
%   Constraint 1:
%
%       r_t = C*z_t 
%   =>  z_t = pinv(C)*r_t + M*s_t,
%
%   where
%
%       M = (I - pinv(C)*C)
%
%   and s_t is some vector T.B.D. by a future constraints.  The first term
%   in Constraint 1 is determined by r_t through the rangespace of C; the
%   second term is determined by s_t via the nullspace of C (notice that 
%   C*M = 0).
%
%   Constraint 2:
%
%       r_t = C*z_t
%           = C*A*z_{t-1} + C*B*u_{t-1} 
%           = C*A*z_{t-1} + 0
%           = C*A*(pinv(C)*r_{t-1} + M*s_{t-1})
%   => s_{t-1} = pinv(C*A*M)*(r_t - C*A*pinv(C)*r_{t-1})
%
%   Constraint 3: for t > 1
%   
%       r_t = C*z_t
%           = C*A*z_{t-1} + C*B*u_{t-1}
%           = C*A*z_{t-1} + 0
%           = C*A*(A*z_{t-2} + B*u_{t-2})
%           = C*A^2*z_{t-2} + C*A*B*u_{t-2}
%   =>  u_{t-2} = pinv(C*A*B)*(r_t - C*A^2*z_{t-2})
%
%   Then after some mindless substitution of Constraints 1 and 2 into 3, we
%   get:
%
%       zstar_t = P1*r_t + P2*r_{t+1}
%
%       ustar_t = Q3*r_{t+2} + Q1*zstar_t
%
%   for the matrices P1,P2,Q1,Q3 defined below.  This allows us to compute
%   the controls for all time without a loop.


% check relative degree
L = 0; q = 0;
AA = eye(size(A));

while all(gather(L==0))
    L = C*AA*B;
    AA = A*AA;
    q=q+1;
end
if ~(q==2)
    error('the relative degree of this system is not equal to 2 -- jgm');
end

% define
Ncases = size(R,1);
Ninputs = size(B,2);
T = size(R,3);
Nstates = size(C,2);
M = eye(Nstates) - pinv(C)*C;



% compute desired states for all time ("one fell swoop")
P1 = (eye(Nstates) - M*pinv(C*A*M)*C*A)*pinv(C);
P2 = M*pinv(C*A*M);
Zstar = shortdata(Ncases,3,(P1*longdata(R(:,:,1:T-1))')') +...
    shortdata(Ncases,3,(P2*longdata(R(:,:,2:T))')');

% compute the corresponding control for all time
Q1 = -pinv(L)*C*A^2;
Q3 = pinv(L);
Ustar = zeros([Ncases,Ninputs,T],'like',R);
Ustar(:,:,1:T-2) = shortdata(Ncases,3,(Q3*longdata(R(:,:,3:T))')') +...
    shortdata(Ncases,3,(Q1*longdata(Zstar(:,:,1:T-2))')');


% the last Zstar and the last two U's are really undetermined, but we just
% filled in the controls with zeros.  For the desired state we just repeat
% the last one.
Zstar(:,:,T) = Zstar(:,:,T-1);



end
%-------------------------------------------------------------------------%






