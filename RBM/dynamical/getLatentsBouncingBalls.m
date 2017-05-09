function [X,Q] = getLatentsBouncingBalls(Nexamples,yrclass,T,balls,m,varargin)
% getLatentsBouncingBalls
%
% USAGE:
%   [X,Q] = getLatentsBouncingBalls(Nexamples,yrclass,T,params.balls,...
%       ones(1,params.balls.N),params.N)

%-------------------------------------------------------------------------%
% Revised: 01/28/16
%   -vectorized, speeding up by an order of magnitude
%   -corrected bug, changing norm used in velocities from matlab's
%   "(matrix) two-norm" to the Frobenius norm, sqrt(sum(v(:).^2)).
% Revised: 12/31/16 
% Revised: ??/??/16
%   -vectorized, rationalized
% Cribbed: ??/??/16
%   from Ilya Sutskever's code
%-------------------------------------------------------------------------%

% params
T = defaulter('sequencelength',T,varargin{:});
Ntraj = floor(Nexamples/T);
Nballs = balls.N;
dt = 0.5;

% malloc
Z = zeros(Ntraj,2,Nballs,T);
V = zeros(Ntraj,2,Nballs,T);
z = zeros(Ntraj,2,Nballs);

% initialize velocity
v = randn(Ntraj,2,Nballs);
v = balls.speed*(v./sqrt(sum(v(:,:).^2,2)));
oo = ones(1,balls.N);

for iTraj = 1:Ntraj

    % initialize position
    good_config=false;
    while ~good_config        
        z(iTraj,:,:) = scalefxn(rand(2,balls.N),0*oo,1*oo,balls.r,...
            balls.boundingbox-balls.r);
        good_config=true;
        
        for iBall=1:balls.N
            for jBall=1:(iBall-1)
                if norm(z(iTraj,:,iBall) - z(iTraj,:,jBall)) <...
                        (balls.r(iBall) + balls.r(jBall))
                    good_config=false;
                end
            end
        end
    end
end

r = permute(balls.r,[3,1,2]);
for t=1:T
    % for how long do we show small simulation
    Z(:,:,:,t) = z;
    V(:,:,:,t) = v;
    
    %%% Wtf?  Don't show every frame, ok; but why show frames inversely as 
    %%% often as they are simulated??
    for mu = 1:round(1/dt)
        
        z = z + dt*v;
        badinds = (z-r) < 0;
        v(badinds) = abs(v(badinds));
        badinds = (z+r) > balls.boundingbox;
        v(badinds) = -abs(v(badinds));
        
        for iBall=2:balls.N
            for jBall=1:(iBall-1)
                
                % how far are these balls from each other?
                ballvecs = diff(z(:,:,[jBall,iBall]),[],3);
                balldsts = sqrt(sum(ballvecs.^2,2));
                
                % the bouncing off part
                w = ballvecs./balldsts;
                v_i = sum(w.*v(:,:,iBall),2);
                v_j = sum(w.*v(:,:,jBall),2);
                [new_v_i, new_v_j] = new_speeds(m(iBall), m(jBall), v_i, v_j);
                velinc_i = w.*(new_v_i - v_i);
                velinc_j = w.*(new_v_j - v_j);
                
                % update the balls that are actually too close
                badTrajs = balldsts < balls.r(iBall)+balls.r(jBall);
                v(badTrajs,:,iBall) = v(badTrajs,:,iBall) + velinc_i(badTrajs,:);
                v(badTrajs,:,jBall) = v(badTrajs,:,jBall) + velinc_j(badTrajs,:);
                
            end
        end
    end
end

if strcmp(yrclass,'gpuArray')
    Z = gpuArray(Z);
    V = gpuArray(V);
end

X           = longdata(reshape(cat(2,Z,V),[Ntraj,4*Nballs,T]));
Q.restarts  = 1:T:(Ntraj*T);
Q.T         = T;

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [new_v1,new_v2] = new_speeds(m1, m2, v1, v2)

new_v2 = (2*m1*v1 + v2*(m2-m1))/(m1+m2);
new_v1 = new_v2 + (v2 - v1);

end
%-------------------------------------------------------------------------%