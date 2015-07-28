function LDSparams = EM4LDS(Nstates,params,varargin)
% Expectation-maximization for the Kalman filter
%
% USAGE: 
%   LDSparams = EM4KF(Nstates,T,params);
%   LDSparams = EM4KF(Nstates,T,params,LDSparams);
%
% This function learns the Kalman filter parameters using
% expectation-maximization.  Currently, it runs in "batch" fashion, using
% Ncases copies of the trajectories to train on before updating the
% parameters.

%-------------------------------------------------------------------------%
% Revised: 08/26/14
%   -eliminated "noisily observed controls"
% Revised: 04/01/14
%   -forced the inference algorithms to use prestored SigmaV & SigmaY---if 
%   they exist
%   -added plotting function for log-likelihoods
%   -added log likelihood (negative cross entropy) for p(v)
% Revised: 03/10/14
%   -made some adjustments for 1D data
% Revised: 01/10/14
%   -rewrote initLDSparams to accommodate controls; fixed some of its bugs
% Revised: 01/07/14
%   -eliminated outputs eKF and eRTSS
% Revised: 12/02/13
%   -changed EM loop to hold off update until after all Ncases
%   trajectories; then added outer loop through "epochs"
% Revised: 11/25/13
%   -adapted gatherXpctSuffStats to allow for nonzero residual means in the
%   emission (something like the column of ones...)
% Revised: 10/30/13
% Created: 10/21/13
%   by JGM
%-------------------------------------------------------------------------%


% Ns
Ncases = params.Ncases;
Nepochs = params.DBNmaxepoch;
TOPLOT = 0;

% malloc
aLL = zeros(1,Nepochs);

% prepare figure
if TOPLOT
    fignum = 456;
    [subplotHandleY,dataplotHandleY] = prepareFigure(fignum,2,3);
    % ppp = rank(LDSparams.C);
    %%% just look at as many states as you can pick up from Y
    % [subplotHandleX,dataplotHandleX] = prepareFigure(2,ppp,3);
    figureInfo.num = fignum;
    figureInfo.subplotHandle = subplotHandleY;
    figureInfo.dataplotHandle = dataplotHandleY;
else
    figureInfo = [];
end
[subplotHandleLL,dataplotHandleLL] = prepareFigure(34,1,2);


% loop through epochs....
for iEpoch = 1:Nepochs
    
    % generate fresh data
    if mod(iEpoch,5)==1
        LDSdata = getLDSdata(params); 
        Y = LDSdata.Y;
        try U = shiftdim(LDSdata.U,1); 
        catch ME, U = zeros(1,size(Y,3),Ncases); end
        try SigmaY = shiftdim(LDSdata.SigmaY,1);        % Time-varying?
        catch ME, SigmaY = LDSparams.SigmaY; end
    end
    
    % initialize model params
    if iEpoch==1
        if isempty(varargin)
            LDSparams = initLDSparams(LDSdata,Nstates,params.dynamics);
        else
            LDSparams = varargin{1};
        end
        aLLold = -Inf;
    end
    
    %------------------------------- E STEP ------------------------------%
    [xpctStats,LLy] = Estep4LDS(LDSparams,shiftdim(Y,1),U,SigmaY,...
        figureInfo);
    xpctT = gatherXpctSuffStats(xpctStats,Y);   % more minimal
    if strcmp(params.dynamics.meta, 'RandInitWithEC')
        xpctT = enforceNoisilyObservedControls(xpctT,params.dynamics);
    end
    aLL(iEpoch) = LLy/Ncases;                           % avg. log-lkhd
    %---------------------------------------------------------------------%
    
    % quit if you've converged in average log likelihood
    aLLnew = aLL(iEpoch);
    fractionalImprovement = (aLLnew - aLLold)/abs(aLLold);
    aLLold = aLLnew;
    fprintf('average LL = %f on EM epoch %i\n',aLLnew,iEpoch);
    animatePlot(34,subplotHandleLL,dataplotHandleLL,2:iEpoch,aLL(2:iEpoch));
    
    thr = 0.0005;	% 1/20 of a percent 
    if (fractionalImprovement < thr) && (mod(iEpoch,5)~=1)
        fprintf('exiting with delta LL below tolerance\n');
        break
    end
    
    
    %----------------------------- M STEP --------------------------------%
    LDSparams = getMLparams(xpctT,LDSparams);
    %---------------------------------------------------------------------%
    
    
end

LDSparams.LL = aLLnew;
fprintf('average LL = %f on EM epoch %i\n',aLLnew,iEpoch);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function LDSparams = initLDSparams(LDSdata,Nx,dynamics)


% unload
Y = LDSdata.Y;      
Ny = size(Y,2);
if isfield(LDSdata,'U'), U = LDSdata.U; Nu = size(U,2); end

% Ns
[Ncases,Nobsvs,T] = size(Y);

switch dynamics.meta
    
    case 'InitWithData'
        
        % data
        Y0 = Y(:,:,1)';
        Ybroad = longdata(Y)';
        
        % initialize C using PCA (assumes no output noise)
        C = estimateOutputMatrix(Nx,Ybroad);
        LDSparams.C = C;
        
        % for building A, you'll need data that have as much rank as Nx
        xfakeshort = shortdata(Ncases,3,(pinv(C)*Ybroad)');
        xfakeshort = xfakeshort(:,1:rank(C),:);
        Xd = xfakeshort;
        while rank(longdata(xfakeshort)) < Nx
            Xd = diff(Xd,[],3);
            xfakeshort = cat(2,xfakeshort(:,:,1:end-1),Xd);
        end
        
        % init A based on data inferred via (initial) C
        xfakef = longdata(xfakeshort(:,:,2:end))';
        if isfield(LDSdata,'U')
            xfakep = [longdata(xfakeshort(:,:,1:end-1))';...
                U(:,:,1:end-1);...
                ones(1,(size(xfakeshort,3)-1)*Ncases,'like',Y)];
            betaTRANS = (xfakef*xfakep')/(xfakep*xfakep');
            LDSparams.A = betaTRANS(:,1:Nx) +...
                randn(Nx)/50;
            LDSparams.B = betaTRANS(:,(Nx+1):(Nx+Nu)) +...
                randn(Nx,Nu)/50;
        else
            xfakep = [longdata(xfakeshort(:,:,1:end-1))';...
                ones(1,(size(xfakeshort,3)-1)*Ncases,'like',Y)];
            betaTRANS = (xfakef*xfakep')/(xfakep*xfakep');
            LDSparams.A = betaTRANS(:,1:end-1);
        end
        LDSparams.muX = betaTRANS(:,end);
        LDSparams.SigmaX = (xfakef*xfakef') - betaTRANS*(xfakep*xfakef');
        
        % init ICs based on initial observations
        foo = randn(Nx,Nx,'like',Y)/1000;         % small & random
        LDSparams.Info0 = pinv(cov((pinv(C)*Y0)')) + foo'*foo; % so it's invertible!
        LDSparams.mu0 = mean(pinv(C)*Y0,2);
        
        % tell getMLparams that muX may be nonzero (see NOTE below)
        % LDSparams.muX = zeros(Nstates,1,'like',Y);
        % LDSparams.muY = zeros(Ndims,1,'like',Y);
        
        
        
    case 'RandInit'

        LDSparams.A = randn(Nx,Nx,'like',Y);
        if isfield(LDSdata,'U'), LDSparams.B = randn(Nx,Nu,'like',Y); end
        foo = randn(Nx,Nx,'like',Y);
        LDSparams.SigmaX = foo'*foo; % 1.0e-05*[0.5 0; 0 0.5];
        LDSparams.muX = randn(Nx,1,'like',Y); % [0;0];
        LDSparams.C = randn(Ny,Nx,'like',Y); % [1 0];
        foo = randn(Nx,Nx,'like',Y);
        LDSparams.Info0 = foo'*foo; % inv([0.0030 0; 0 0.1000]);
        LDSparams.mu0 = randn(Nx,1,'like',Y); % [1.0236; 2.5];
        
        
    case 'RandInitWithEC'
        
        [Ny,Nx] = size(dynamics.C);
        [Nv,Nec] = size(dynamics.H);
        
        muU = randn(Nec,1,'like',Y);
        foo = randn(Nec,Nec,'like',Y);
        InfoU = foo'*foo;
        
        LDSparams.A=[randn(Nx,Nx+Nec,'like',Y);zeros(Nec,Nx+Nec,'like',Y)];
        if isfield(LDSdata,'U')
            LDSparams.B = [randn(Nx,Nu,'like',Y); zeros(Nec,Nu,'like',Y)];
        end
        LDSparams.muX = [randn(Nx,1,'like',Y); muU];
        foo = randn(Nx,Nx,'like',Y);
        LDSparams.SigmaX = [foo'*foo, zeros(Nx,Nec,'like',Y);...
            zeros(Nec,Nx,'like',Y), inv(InfoU)];
        LDSparams.C = [randn(Ny,Nx,'like',Y), zeros(Ny,Nec,'like',Y);....
            zeros(Nv,Nx,'like',Y),randn(Nv,Nec,'like',Y)];
        LDSparams.mu0 = [randn(Nx,1,'like',Y); muU];
        foo = randn(Nx,Nx,'like',Y);
        LDSparams.Info0 = inv([inv(foo'*foo), zeros(Nx,Nec,'like',Y);...
            zeros(Nec,Nx) inv(InfoU)]);
        LDSparams.T = T;
        
        
end
   

% make the output and transition noises small [wrt what??]
% foo = randn(Ndims,'like',Y)/1000;
% LDSparams.SigmaY = foo'*foo;
%%% alternatively, for systems with nonconstant covariance, you could just
%%% give this param to the filter:
%%% LDSparams.SigmaY(1:Ndims,1:Ndims,1:T) = squeeze(LDSdata.SigmaY(1,:,:,:,:));
LDSparams.T = T;


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function P = estimateOutputMatrix(Ncols,Y)
% Y = P*X, where size(P) = Nrows x Ncols.  Estimate P from the intrinsic
% (linear) dimensionality of Y.

% params
Nrows = size(Y,1);
NpcsMAX = 10;
%%% you could really use a better criterion for max num pcs than just 10.

% Npcs <= Nrows, Npcs <= Ncols, Npcs <= NpcsMAX
Npcs = min([NpcsMAX,Nrows,Ncols]);  

% how many meaningful dimensions are there in Y?
[V,D] = eig(cov(Y'));
[maxima,indices] = sort(diag(D),'descend');     % sort by eigenvalue
W = V(:,indices(1:Npcs));                       % select 1st Npcs PCs
P = zeros(Nrows,Ncols,'like',Y);                % a P that does nothing
P(:,1:Npcs) = W;                                % the useful part of P

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [subplotHandle,dataplotHandle] = prepareFigure(fignum,Ndims,Ndatavecs)

figure(fignum); clf; hold on;
for dim = 1:Ndims
    subplotHandle = subplot(1,Ndims,dim);
    hold on;
    for j = 1:Ndatavecs
        dataplotHandle(dim,j) = plot(NaN,NaN);
    end
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [xpctStats,LLy] = Estep4LDS(LDSparams,Y,U,SigmaY,figureInfo)


% filter and smooth
KF = arrayfun(@(iCase)(KalmanFilter(setfield(LDSparams,'SigmaY',...
    SigmaY(:,:,:,iCase)),Y(:,:,iCase),U(:,:,iCase))),1:size(Y,3));
RTSS = arrayfun(@(iCase)(RTSsmoother(setfield(LDSparams,...
    'SigmaY',SigmaY(:,:,:,iCase)),KF(iCase))),1:size(Y,3));

% store results
xpctStats.xpctX = shiftdim(cat(3,RTSS(:).XHAT),2);
xpctStats.AvgCvrnX = mean(cat(4,RTSS(:).XCVRN),4);
xpctStats.AvgCvrnXfXp = mean(cat(4,RTSS(:).XfXpCVRN),4);
LLy = -sum([KF(:).XNtrpY]);


% plot?
if ~isempty(figureInfo)
    iCase = 5;
    % for iCase = 1:size(Y,3)
        % the only thing you can really test: outputs!!
        shatKF = LDSparams.C*KF(iCase).XHATMU;
        shatRTSS = LDSparams.C*RTSS(iCase).XHAT;
        y = Y(:,:,iCase);
        
        % now write to the plot
        %%%%% hardcoded fignum = 1
        animatePlot(figureInfo.num,figureInfo.subplotHandle,...
            figureInfo.dataplotHandle,1:size(y,2),y,shatKF,shatRTSS);
        % animatePlot(2,subplotHandleX,dataplotHandleX,1:T,zeros(ppp,T),... % x(1:ppp,:),...
        %     XhatKF(1:ppp,:),RTSSdstrbs.XHAT(1:ppp,:));
        % keyboard
        pause();
    % end
end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function animatePlot(fignum,subplotHandle,dataplotHandle,xdata,varargin)

% params
Ndatavecs = length(varargin);

% pull up the current figure
set(0,'CurrentFigure',figure(fignum))
colors = 'krbcmy';      
%%% try not to use any more than that

for ivec = 1:Ndatavecs
    
    ydata = varargin{ivec};
    Ndims = size(ydata,1);
    
    for dim = 1:Ndims
        set(dataplotHandle(dim,ivec),'XData',xdata,'YData',ydata(dim,:),...
            'color',colors(ivec));
    end
end
%%% set(subplotHandle,'XLim',[0 Max]);
%%% if you want to keep the x or y limits fixed over the animation....

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function xpctT = gatherXpctSuffStats(xpctStats,Y)
% NB that the sufficient statistics are all in the form of *average*
% (across time and cases) expected outer products, and of course that:
%
%   <E[X*X'|y]> = <E[X|y]E[X'|y]> + <Cov[X|y]>
%
%   but
%
%   <E[X*y'|y]> = <E[X|y]Y'>


%                   -------------NOTE--------------
% There's a redundancy here.  A translation of Y (muY!=0) can usually be
% absorbed into a constant input to X, namely muX := (I-A)C'inv(CC')muY.
% Proceeding backwards from the conclusion:
%
%       Z_{t+1} = A*Z_t + (I-A)*C'*inv(C*C')*muY
%       Y_t = C*Z_t
%   =>  Z_{t+1} = A*(Z_t - C'*inv(C*C')*muY) + C'*inv(C*C')*muY
%   =>  Z_{t+1} - C'*inv(C*C')*muY = A*(Z_t - C'*inv(C*C')*muY)
%   X := Z - C'*inv(C*C')*muY => 
%       X_{t+1} = A*X_t.
%       Y_t = C*(X_t + C'*inv(C*C')*muY)
%           = C*X_t + muY.
%
% NB: *But this assumes that C*C' is invertible* ("C is fat").
%
% Likewise, a constant input to X (muX!=0) can be absorbed into a 
% translation of the emission Y, namely muY = C*inv(I-A)*muX:
% 
%       X_{t+1} = A*X_t
%       Y_t = C*X_t + C*inv(I-A)*muX
%           = C*(X_t + inv(I-A)*muX)
%   Z := X + inv(I-A)*muX => 
%           = C*Z_t
%   =>  Z_{t+1} = X_{t+1} + inv(I-A)*muX
%               = A*X_t + inv(I-A)*muX
%               = A*(Z_t - inv(I-A)*muX) + inv(I-A)*muX
%               = A*Z_t + (I-A)*inv(I-A)*muX
%               = A*Z_t + muX.
%
% But this assumes that (I-A) is invertible, which is much less likely!!  
% (It is not the case for your default dynamics, for example).  (You might
% think you could pull a similar trick to the previous by using the
% pseudoinverse of (I-A), but rank(Q*Q') = rank(Q) for real matrices.)
% 
% Therefore, by default this function will assume that the only bias is on
% the transitions.  ***HOWEVER***, you should check before using this that
%
% 	rank(C*C') = Nstateobsv (won't be true, e.g., for neural array data)
%
% When these are not satisfied you should allow offsets biases for Y.
%                   ------------------------------


% function
mucov2op = @(muA,muB,AvgSgm)(muA*muB'/size(muA,2) + AvgSgm);

% FLAGS
BIASEDTRANSITIONS = 1;  % isfield(LDSparams,'muX');
BIASEDEMISSIONS = 0;    % isfield(LDSparams,'muY');


% sift out data for regression
muX0 = xpctStats.xpctX(:,:,1)';
muXp = longdata(xpctStats.xpctX(:,:,1:end-1))';
muXf = longdata(xpctStats.xpctX(:,:,2:end))';
muX = longdata(xpctStats.xpctX)';
Ybroad = longdata(Y)';

AvgCvrnX = mean(xpctStats.AvgCvrnX,3);
AvgCvrnX0 = xpctStats.AvgCvrnX(:,:,1);
AvgCvrnXpXp = mean(xpctStats.AvgCvrnX(:,:,1:end-1),3);
AvgCvrnXfXf = mean(xpctStats.AvgCvrnX(:,:,2:end),3);
AvgCvrnXfXp = mean(xpctStats.AvgCvrnXfXp,3);

%%%% if there's a (visible) control....

% STATE/STATE AVERAGE EXPECTED OUTER PRODUCTS
% ...over the first state
xpctT.mu0 = mean(muX0,2);
xpctT.x0x0 = mucov2op(muX0,muX0,AvgCvrnX0);

% ...over all but the first, all but the last, and their interaction
xpctT.XpXp = mucov2op(muXp,muXp,AvgCvrnXpXp);
xpctT.XfXf = mucov2op(muXf,muXf,AvgCvrnXfXf);
xpctT.XfXp = mucov2op(muXf,muXp,AvgCvrnXfXp);
if BIASEDTRANSITIONS
    xpctT.XpXp = [xpctT.XpXp mean(muXp,2); mean(muXp,2)' 1];
    xpctT.XfXp = [xpctT.XfXp mean(muXf,2)];
end

% STATE/EMISSION AVERAGE EXPECTED OUTER PRODUCTS
xpctT.XX = mucov2op(muX,muX,AvgCvrnX);
xpctT.YX = mucov2op(Ybroad,muX,0);
xpctT.YY = mucov2op(Ybroad,Ybroad,0);
if BIASEDEMISSIONS
    xpctT.XX = [xpctT.XX mean(muX,2); mean(muX,2)' 1];
    xpctT.YX = [xpctT.YX mean(Ybroad,2)];
end


end
%-------------------------------------------------------------------------%