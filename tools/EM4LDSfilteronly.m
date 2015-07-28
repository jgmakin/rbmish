function LDSparams = EM4LDSfilteronly(Nstates,Ninputs,T,params,varargin)
% Expectation-maximization for the Kalman filter
%
% USAGE: LDSparams = EM4KF(T,params,learningtype)
%
% This function learns the Kalman filter parameters using
% expectation-maximization.  Currently, it runs in "batch" fashion, using
% Ncases copies of the trajectories to train on before updating the
% parameters.


%%% It doesn't *really* make sense, in general, to plot the hidden states,
%%% since you're really just guessing that the first two states are
%%% important.


%-------------------------------------------------------------------------%
% Revised: 04/01/14
%   -forced the inference algorithms to use prestored SigmaV and
%   SigmaY---if they exist
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
Ndims = params.Ndims;
Ncases = params.Ncases;
maxepochs = 90;
TOPLOT = 0;

% prepare figure
if TOPLOT
    [subplotHandleY,dataplotHandleY] = prepareFigure(1,Ndims,3);
    % ppp = rank(LDSparams.C);
    %%% just look at as many states as you can pick up from Y
    % [subplotHandleX,dataplotHandleX] = prepareFigure(2,ppp,3);
end
[subplotHandleLL,dataplotHandleLL] = prepareFigure(34,1,2);

% malloc
xpctX = zeros(Ncases,Nstates,T);
xpctU = zeros(Ncases,Ninputs,T);
aLLy = zeros(1,maxepochs);
aLLv = zeros(1,maxepochs);


% loop through epochs....
for epoch = 1:maxepochs
    
    % generate fresh data
    if mod(epoch,5)==1
        LDSdata = getLDSdata(T,params);
        Y = LDSdata.Y;
        V = LDSdata.V;
    end
    
    % initialize model params
    if epoch==1
        if isempty(varargin)
            LDSparams = initLDSparams(Nstates,Ninputs,LDSdata);
        else
            LDSparams = varargin{1};
        end
        aLLold = -Inf;
    end
    
    % init cumulative sums
    CUMcvrnX = 0; CUMCvrnXfXp = 0; CUMcvrnU = 0; CUMCvrnXfUp = 0;
    CUMCvrnXpUp = 0; LLy = 0; LLv = 0;
    
    %------------------------------- E STEP ------------------------------%
    % loop through trajectories
    for iCase = 1:Ncases
        
        % emissions for this trajectory
        y(1:Ndims,1:T) = squeeze(Y(iCase,:,:,:));
        if isfield(LDSdata,'SigmaY')
            %%% just tell it what the emission variances are
            LDSparams.SigmaY(1:Ndims,1:Ndims,1:T) =...
                squeeze(LDSdata.SigmaY(iCase,:,:,1,:));
        end
        if Ninputs>0
            LDSparams.V(1:Ndims,1:T) = squeeze(V(iCase,:,:,:)); 
            if isfield(LDSdata,'SigmaV')
                %%% just tell it what the emission variances are
                LDSparams.SigmaV(1:Ndims,1:Ndims,1:T) =...
                    squeeze(LDSdata.SigmaV(iCase,:,:,1,:));
            end
        end
        %%% at the moment, none of the "adjusters" is employed (see
        %%% KF4PPC.m), which are required for bouncing, etc.
        
        % run filter, smoother
        KFdstrbs = KalmanFilter(LDSparams,y);
        %%%%%%%
        for t = 1:(T-1)
            Pt = KFdstrbs.CVRNMU(:,:,t);
            Jt = Pt*LDSparams.A'*KFdstrbs.INFOTU(:,:,t+1);
            XfXpCVRN(:,:,t) = KFdstrbs.CVRNMU(:,:,t+1)*Jt';            
            RTSSdstrbs.XfXpCVRN(:,:,t) =...
                Pt + Jt*XfXpCVRN(:,:,t) - Jt*LDSparams.A*Pt;
        end        
        RTSSdstrbs.XHAT = KFdstrbs.XHATMU;
        RTSSdstrbs.XCVRN = KFdstrbs.CVRNMU;
        %%%%%%%
            
        % plot
        if TOPLOT
            % the only thing you can really test: outputs!!
            shatKF = LDSparams.C*KFdstrbs.XHATMU;
            shatRTSS = LDSparams.C*RTSSdstrbs.XHAT;
            
            % now write to the plot
            animatePlot(1,subplotHandleY,dataplotHandleY,1:T,y,shatKF,shatRTSS);
            % animatePlot(2,subplotHandleX,dataplotHandleX,1:T,zeros(ppp,T),... % x(1:ppp,:),...
            %     XhatKF(1:ppp,:),RTSSdstrbs.XHAT(1:ppp,:));
            pause();
        end
        
        % accumulate
        xpctX(iCase,:,:) = RTSSdstrbs.XHAT;
        CUMcvrnX = CUMcvrnX + RTSSdstrbs.XCVRN;
        CUMCvrnXfXp = CUMCvrnXfXp + RTSSdstrbs.XfXpCVRN;
        LLy = LLy - KFdstrbs.XNtrpY;
        if Ninputs>0
            fprintf('.');
            xpctU(iCase,:,:) = RTSSdstrbs.UHAT;
            CUMcvrnU = CUMcvrnU + RTSSdstrbs.UCVRN;
            CUMCvrnXfUp = CUMCvrnXfUp + RTSSdstrbs.XfUpCVRN;
            CUMCvrnXpUp = CUMCvrnXpUp + RTSSdstrbs.XpUpCVRN;
            LLv = LLv - KFdstrbs.XNtrpV;
        end
    end
    aLLy(epoch) = LLy/Ncases;
    aLLv(epoch) = LLv/Ncases;
    
    
    % gather expected statistics
    xpctStats.xpctX = xpctX;
    xpctStats.AvgCvrnX = CUMcvrnX/Ncases;
    xpctStats.AvgCvrnXfXp = CUMCvrnXfXp/Ncases;
    if Ninputs>0
        xpctStats.xpctU = xpctU;
        xpctStats.AvgCvrnU = CUMcvrnU/Ncases;
        xpctStats.AvgCvrnXfUp = CUMCvrnXfUp/Ncases;
        xpctStats.AvgCvrnXpUp = CUMCvrnXpUp/Ncases;
    end
    
    % turn into (more) minimal sufficient stats
    xpctT = gatherXpctSuffStats(xpctStats,Y,V);
    %---------------------------------------------------------------------%
    
    
    % quit if you've converged in average log likelihood
    aLL = aLLy(epoch); % + aLLv;
    fractionalImprovement = (aLL - aLLold)/abs(aLLold);
    aLLold = aLL;
    fprintf('average LLy = %f, LLv = %f on EM epoch %i\n',...
        aLLy(epoch),aLLv(epoch),epoch);
    animatePlot(34,subplotHandleLL,dataplotHandleLL,1:epoch,...
        aLLy(1:epoch),aLLv(1:epoch));
    
    thr = 0.0005; % 0.005;   % half a percent
    if fractionalImprovement < thr   
        fprintf('exiting with delta LL below tolerance\n');
        break
    end
    
    
    %----------------------------- M STEP --------------------------------%
    LDSparams = getMLparams(xpctT,LDSparams);
    %---------------------------------------------------------------------%
    
    
end

LDSparams.LL = aLL;
fprintf('average LL = %f on EM epoch %i\n',aLL,epoch);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function LDSparams = initLDSparams(Nstates,Ninputs,LDSdata)


% unload structure
Y = LDSdata.Y;

% Ns
[Ncases,Ndims,Nmods,T] = size(Y);
Y0 = Y(:,:,:,1)';
Ybroad = longdata(Y)';

% initialize C using PCA (assumes no output noise)
C = estimateOutputMatrix(Nstates,Ybroad);
LDSparams.C = C;
LDSparams.C = randn(Ndims,Nstates); % [1 0]; overwrite with randomness

% for controlled systems, init H, SigmaV, muV, and prior over U
if Ninputs>0
    
    Vbroad = longdata(LDSdata.V)';
    H = estimateOutputMatrix(Ninputs,Vbroad);
    LDSparams.H = H;
    LDSparams.H = randn(Ninputs,Ndims); % 1; overwrite with randomness
       
    % make the control-observation noise small [wrt what??]
    % foo = randn(Ndims);                         % random (/100?)
    % LDSparams.SigmaV = foo'*foo;
    %%% alternatively, for systems with nonconstant covariance, you could
    %%% just give these params to the filter: 
    %%% LDSparams.SigmaV(1:Ndims,1:Ndims,1:T) = squeeze(LDSdata.SigmaV(1,:,:,:,:));

        
    % intialize prior over U based on (initial) H
    foo = randn(Ninputs)/100;                   % small & random
    LDSparams.InfoU = pinv(cov((pinv(H)*Vbroad)')) + foo'*foo; % invertible!
    LDSparams.muU = mean(pinv(H)*Vbroad,2);    
    foo = randn(Ninputs);                       % overwrite with randomness
    LDSparams.InfoU = foo'*foo/1000; % 3.656e2;
    LDSparams.muU =  randn(Ninputs,1); % -0.001754435768020;
end


% for building A, you'll need data that have as much rank as Nstates
xfakeshort = shortdata(Ncases,4,(pinv(C)*Ybroad)');
xfakeshort = xfakeshort(:,1:rank(C),:,:);
Xd = xfakeshort;
while rank(longdata(xfakeshort)) < Nstates
    % Xd = Xd(:,:,:,2:end);
    Xd = diff(Xd,[],4);
    xfakeshort = cat(2,xfakeshort(:,:,:,1:end-1),Xd);
end
% xfakeshort = shortdata(Ncases,3,(pinv(C)*Ybroad)');

% init A based on data inferred via (initial) C (and H if controlled)
xfakef = longdata(xfakeshort(:,:,:,2:end))';
if Ninputs > 0
    ufakep = pinv(H)*longdata(LDSdata.V(:,:,:,1:(size(xfakeshort,4)-1)))';
    xfakep = [longdata(xfakeshort(:,:,:,1:end-1))'; ufakep;...
        ones(1,(size(xfakeshort,4)-1)*size(xfakeshort,1))];
    
    betaTRANS = (xfakef*xfakep')/(xfakep*xfakep');
    LDSparams.A = betaTRANS(:,1:Nstates) + randn(Nstates)/50;
    LDSparams.B = betaTRANS(:,(Nstates+1):(Nstates+Ninputs)) + randn(Nstates,Ninputs)/50;
    LDSparams.A = randn(Nstates,Nstates);       % overwrite with randomness
    LDSparams.B = randn(Nstates,Ninputs);       % [0; 1];
else
    xfakep = [longdata(xfakeshort(:,:,:,1:end-1))';...
        ones(1,(size(xfakeshort,4)-1)*size(xfakeshort,1))];
    betaTRANS = (xfakef*xfakep')/(xfakep*xfakep');
    LDSparams.A = betaTRANS(:,1:end-1);
    LDSparams.A = randn(Nstates,Nstates);       % overwrite with randomness
end
LDSparams.muX = betaTRANS(:,end);
LDSparams.SigmaX = (xfakef*xfakef') - betaTRANS*(xfakep*xfakef');
foo = randn(Nstates,Nstates);
LDSparams.SigmaX = foo'*foo; % 1.0e-05*[0.5 0; 0 0.5];
LDSparams.muX = randn(Nstates,1); % [0;0];      % overwrite with randomness

% tell getMLparams that muX may be nonzero (see NOTE below)
% LDSparams.muX = zeros(Nstates,1);
% LDSparams.muY = zeros(Ndims,1);

% intialize state 0 based on (initial) C
foo = randn(Nstates)/1000;                      % small & random
LDSparams.Info0 = pinv(cov((pinv(C)*Y0)')) + foo'*foo; % so it's invertible!
LDSparams.mu0 = mean(pinv(C)*Y0,2);
foo = randn(Nstates,Nstates);
LDSparams.Info0 = foo'*foo; % inv([0.0030 0; 0 0.1000]);
LDSparams.mu0 = randn(Nstates,1); % [1.0236; 2.5];overwrite with randomness


% make the output and transition noises small [wrt what??]
% foo = randn(Ndims)/1000;                          % small & random
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
P = zeros(Nrows,Ncols);                         % a P that does nothing
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
function xpctT = gatherXpctSuffStats(xpctStats,Y,V)
% NB that the sufficient statistics are all in the form of *average*
% (across time and cases) expected outer products, and of course that:
%
%   <E[X*X'|y]> = <E[X|y]E[X'|y]> + <Cov[X|y]>
%
%   but
%
%   <E[X*y'|y]> = <E[X|y]Y'>



%%% xpctX,AvgCvrnX,AvgCvrnXpXf,Y)


%                   -------------NOTE--------------
% There's a redundancy here.  A translation of Y (muY!=0) can always be
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
% Something similar is true of the controls.  If H*H' is invertible, then
% we can write muX = -B*H'*inv(H*H')*muV.  Then:
%
%       V_t = H*U_t
%       Z_{t+1} = A*Z_t + B*U_t - B*H'*inv(H*H')*muV
%               = A*Z_t + B*[U_t - H'*inv(H*H')*muV]
%   UU_t := U_t - H'*inv(H*H')*muV =>
%               = A*Z_t + B*UU_t
%   =>  V_t = H*(UU_t + H'*inv(H*H')*muV)
%           = H*UU_t + H*H'*inv(H*H')*muV
%           = H*UU_t + muV.
%
% For you, H will generally be the identity, so H*H' will be invertible.
% NB: if H is "tall" (the controls, or combinations thereof, appear
% redundantly in V), then this will not work, and you need to allow a bias
% muV!!
%
% Therefore, by default this function will assume that the only bias is on
% the transitions.  ***HOWEVER***, you should check before using this that
%
% (1) rank(C*C') = Nstateobsv (won't be true, e.g., for neural array data)
% (2) rank(H*H') = Nctrlobsvs
%
% When these are not satisfied you should allow offsets biases for Y and V.
%                   ------------------------------



% function
mucov2op = @(muA,muB,AvgSgm)(muA*muB'/size(muA,2) + AvgSgm);

% FLAGS
BIASEDTRANSITIONS = 1;  % isfield(LDSparams,'muX');
BIASEDEMISSIONS = 0;    % isfield(LDSparams,'muY');
BIASEDCTRLOBSVS = 0;    % isfield(LDSparams,'muV');
CTRLD = isfield(xpctStats,'xpctU');



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

if CTRLD
   
    % sift out data
    muUp = longdata(xpctStats.xpctU(:,:,1:end-1))';
    muU = longdata(xpctStats.xpctU)';
    AvgCvrnU = mean(xpctStats.AvgCvrnU,3);
    AvgCvrnUpUp = mean(xpctStats.AvgCvrnU(:,:,1:end-1),3);
    AvgCvrnXpUp = mean(xpctStats.AvgCvrnXpUp,3);
    AvgCvrnXfUp = mean(xpctStats.AvgCvrnXfUp,3);
    Vbroad = longdata(V)';
    
    % CONTROL/CONTROL OBSV. AVERAGE (over time) EXPECTED OUTER PRODUCTS
    xpctT.UU = mucov2op(muU,muU,AvgCvrnU);
    xpctT.VU = mucov2op(Vbroad,muU,0);
    xpctT.VV = mucov2op(Vbroad,Vbroad,0);
    xpctT.muU = mean(muU,2);
    if BIASEDCTRLOBSVS
        xpctT.UU = [xpctT.UU mean(muU,2); mean(muU,2)' 1];
        xpctT.VU = [xpctT.VU mean(Vbroad,2)];
    end
    
    % muXp -> [muXp; muUp]    
    muXp = [muXp; muUp];
    AvgCvrnXpXp = [AvgCvrnXpXp, AvgCvrnXpUp; AvgCvrnXpUp', AvgCvrnUpUp];
    AvgCvrnXfXp = [AvgCvrnXfXp AvgCvrnXfUp];
end



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











