function C = particleIRLS(W,xpctT,LDSparams,dstrb)
% ParticleIRLS  IRLS for particle-based emissions of LTI systems
%
% ParticleIRLS(W,xpctT,LDSparams,dstrb) takes particle weights W; particle
% trajectories Xtu, XtuW, etc., stored in xpctT; the emission parameters C
% (and possibly others) stored in LDSparams; and the type of dstrb; and
% returns an updated estimate of C, using iteratively reweighted least
% squares.
%
% [[size of inputs vars...]]
% 
%
% N.B. For now, this function assumes that muYX is zero

%-------------------------------------------------------------------------%
% Created: 09/09/17
%   by JGM
%-------------------------------------------------------------------------%

% init
thr = 0.0005;
NiterMax = 100;
Xtu = xpctT.Xtu;
XtuW = Xtu.*shiftdim(W,-1);
C = LDSparams.C;
if size(Xtu,1) > size(C,2), C = cat(2,LDSparams.muYX,C); end
Nobsvs = size(C,1);
deltaChat = NaN(size(C),'like',W);
Niter = 0;

tic
fprintf('IRLS loop...\n')
while ~all(abs(gather(deltaChat(:)./C(:)))<thr)&&(Niter<NiterMax)
    
    % (weighted) particles to sufficient statistics
    [negH,grad] = particleSufficientStats(Xtu,XtuW,C,xpctT,LDSparams,dstrb);
    
    % Newton-Raphson: deltaphi = -inv(d^2L/dphi^2)*(dL/dphi)
    for iObsv = 1:Nobsvs
        deltaChat(iObsv,:) = grad(iObsv,:)/negH(:,:,iObsv);
    end
    %%%% since the outputs are independent, it may make more sense to make
    %%%% this the *outer* loop, and the while loop the inner one.
    
    % update
    C = C + deltaChat;
	Niter = Niter+1;
    fprintf('.');
end
fprintf('\n');
toc;




end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [negH,grad] = particleSufficientStats(Xtu,XtuW,C,xpctT,...
    LDSparams,dstrb)
% Gather expected sufficient statistic for IRLS.  Quite complicated, but 
% this avoids all loops
% 
% The master ordering is:
%   (Nstates x Nstates x Nobsvs x Nparticles x T)
%


% Ns
Nobsvs = size(C,1);
[~,Nparticles,T] = size(Xtu);


% which version of IRLS?
switch dstrb
    case 'Poisson'
        %%%%
        % This sometimes generates some Infs (and maybe NaNs)
        %%%%
        
        % E[Y|X]                    (Nobsv x Nparticles x T)
        xpctYgivenX = reshape(exp(C*Xtu(:,:)),[Nobsvs,Nparticles,T]);
        
        % E[Y|X]*X*W                (Nstates x Nobsvs x Nparticles x T)
        xpctYgivenXXW = permute(xpctYgivenX,[4,1,2,3]).*permute(XtuW,[1,4,2,3]);
        
        % gradient:                 (Nobsvs x Nstates)
        %   sum_y{<y*X>_W} - sum_y{<E[y|X]*X'>_W}     
        %   = sum_y{<y*X - E[y|X]*X'>_W}            
        grad = xpctT.Y*permute(sum(XtuW,2),[3,1,2]) -...
            sum(sum(xpctYgivenXXW,4),3)';
        
        % negative Hessian:         (Nstates x Nstates x Nobsvs),
        %   sum_y{<X*D(y)*X'>_W}  
        %   where D = diag(E[Y|X]*X*W)
        negH = sum(sum(...
            permute(xpctYgivenXXW,[5,1,2,3,4]).*...
            permute(Xtu,[1,4,5,2,3]),5),4); 
        % Nstates x Nstates x Nobsvs x Nparticles x T
        
        
    case 'GammaFixedScale'
        % See labnotes for a derivation

        % data
        logY = permute(xpctT.Y,[1,3,2]);
        
        % canonical parameters
        ks = reshape(exp(C*Xtu(:,:)),[Nobsvs,Nparticles,T]);
        logth = log(LDSparams.th);
        
        % E[log(Y)|X], Var[log(Y)|X]
        xpctLogYgivenX = psi(0,ks) + logth;
        vrncLogYgivenX = psi(1,ks);
        
        % (log(Y) - E[log(Y)|X])*k
        v = (logY - xpctLogYgivenX).*ks;
        
        % -Var[log(Y)|X]*k^2  [[ + (log(Y) - E[log(Y)|X])*k ]]
        u = (-vrncLogYgivenX).*(ks.^2);   % + v;
        %%% +v: Hessian vs. expected Hessian
        
        % gradient:                 (Nobsvs x Nstates)
        %   sum_y{<v*X'>_W}
        grad = sum(sum(permute(v,[1,4,2,3]).*permute(XtuW,[4,1,2,3]),4),3);
        
        % negative Hessian:     (Nstates x Nstates x Nobsvs)
        %   sum_y{<X*diag(u(y))*X'>_W}  
        negH = sum(sum(...
            (permute(XtuW,[4,1,5,2,3]).*permute(-u,[4,5,1,2,3])).*...
            permute(Xtu,[1,4,5,2,3]),5),4);
        
        
    otherwise
        error('unexpected emission for particle suff. stats. -- jgm');
end

end
%-------------------------------------------------------------------------%     
