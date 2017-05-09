function LDSparams = EM4LDSfilteronly(Nstates,params,varargin)
% EM *without smoothing* for an LTI system w/Gaussian noise
%
% USAGE:
%   LDSparams = EM4KF(Nstates,params);
%   LDSparams = EM4KF(Nstates,params,'params',LDSparams,'data',LDSdata);
%
% This function learns the Kalman filter parameters using
% expectation-maximization.  Currently, it runs in "batch" fashion, using
% Ncases copies of the trajectories to train on before updating the
% parameters.

%-------------------------------------------------------------------------%
% Revised: 08/26/16
%   -freshly copied over from EM4LDS.m, then edited to again elimate RTS
%   smoother
% Revised: 08/17/15
%   -large revision: rewrote to accommodate LDSdata that has 
%   variable-length trajectories (and hence is an array structure rather
%   than a structure of tensors)
% Revised: 08/26/14
%   -eliminated "noisily observed controls"
% Revised: 04/01/14
%   -forced the inference algorithms to use prestored SigmaV & SigmaYX---if
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
NepochsMax = params.NepochsMax;
DIAGCOVS = 0;
TOPLOT = 0;
Ntest = 5;
thr = 0.0005;	% 1/20 of a percent

% variable arguments
iParams = find(strcmp(varargin,'params'));
iData = find(strcmp(varargin,'data'));
iDiag = find(strcmp(varargin,'diagonal covariances'));
if ~isempty(iDiag), DIAGCOVS = varargin{iDiag+1}; end

% prepare figure
if TOPLOT
    fignum = 456;
    [subplotHandleY,dataplotHandleY] = prepareFigure(fignum,...
        size(params.dynamics.C,1),3);
    % ppp = rank(LDSparams.C);
    %%% just look at as many states as you can pick up from Y
    % [subplotHandleX,dataplotHandleX] = prepareFigure(2,ppp,3);
    figureInfo.num = fignum;
    figureInfo.subplotHandle = subplotHandleY;
    figureInfo.dataplotHandle = dataplotHandleY;
else
    figureInfo = [];
end
[subplotHandleLL,dataplotHandleLL] = prepareFigure(34,1,1);


% loop through epochs....
for iEpoch = 1:NepochsMax
    
    %------------------- DATA GENERATION/INITIALIZATION ------------------%
    % every Ntest epochs...
    if mod(iEpoch,Ntest)==1
        
        if isempty(iData)
            trajs = getLDStrajs(params);
        else
            trajs = varargin{iData+1};
        end
        
        % initialize model params
        if iEpoch==1
            if isempty(iParams)
                LDSparams = initLDSparams(trajs,Nstates,params.dynamics);
            else
                LDSparams = varargin{iParams+1};
            end
            XNtrpYold = -Inf;
            
            % malloc
            XNtrpY = zeros(1,NepochsMax,'like',trajs(1).Y); 
        end
    end
    %---------------------------------------------------------------------%
    
    
    %------------------------------- E STEP ------------------------------%
    [xpctT,XNtrpY(iEpoch)] = Estep4LDS(trajs,LDSparams,figureInfo);
    % if strcmp(pLDSparams.meta, 'RandInitWithEC')
    %%%% broken: you need somehow to pull out separate C and H....
    %     xpctT = enforceNoisilyObservedControls(xpctT,C,H,pLDSparams.T);
    % end
    %---------------------------------------------------------------------%
    
    
    %------------------------ CHECK FOR CONVERGENCE ----------------------%
    % quit if you've converged in average log likelihood
    XNtrpYnew = XNtrpY(iEpoch);
    fractionalImprovement = (XNtrpYold - XNtrpYnew)/abs(XNtrpYold);
    XNtrpYold = XNtrpYnew;
    fprintf('average cross entropy = %f on EM epoch %i\n',XNtrpYnew,iEpoch);
	if usejava('desktop')
    	animatePlot(34,subplotHandleLL,dataplotHandleLL,2:iEpoch,XNtrpY(2:iEpoch));
    end
    if (fractionalImprovement < thr) && (mod(iEpoch,Ntest)~=1)
        fprintf('exiting with delta cross entropy below tolerance\n');
        break
    end
    %---------------------------------------------------------------------%
    
    
    %----------------------------- M STEP --------------------------------%
    LDSparams = ML4LDS(xpctT,LDSparams,'diagonal covariances',DIAGCOVS);
    %---------------------------------------------------------------------%
    
    
end

LDSparams.LL = -XNtrpYnew;
fprintf('average cross entropy = %f on EM epoch %i\n',XNtrpYnew,iEpoch);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function LDSparams = initLDSparams(trajs,Nx,dynamics)

% unload
Ntrajs = length(trajs);
Y0 = arrayfun(@(iCase)(trajs(iCase).Y(:,1)),1:Ntrajs,'UniformOutput',0);
Y0 = cat(2,Y0{:});
Ybroad = cat(2,trajs(:).Y);
Ny = size(Ybroad,1);
if isfield(trajs,'U'), Ubroad = cat(2,trajs(:).U); Nu = size(Ubroad,1); end
Ts = arrayfun(@(iCase)(size(trajs(iCase).Y,2)),1:Ntrajs);
if all(Ts == Ts(1)), T = Ts(1); else T = 'variable'; end


switch dynamics.meta
    
    %%% this still only works with equal-length trajectories
    case 'withData'
        
        % initialize C using PCA (assumes no output noise)
        C = estimateOutputMatrix(Nx,Ybroad);
        LDSparams.C = C;
        
        % for building A, you'll need data that have as much rank as Nx
        xfakeshort = shortdata(Ntrajs,3,(pinv(C)*Ybroad)');
        xfakeshort = xfakeshort(:,1:rank(C),:);
        Xd = xfakeshort;
        while rank(longdata(xfakeshort)) < Nx %%% <= Nx
            Xd = diff(Xd,[],3);
            xfakeshort = cat(2,xfakeshort(:,:,1:end-1),Xd);
        end
        
        % init A based on data inferred via (initial) C
        xfakef = longdata(xfakeshort(:,:,2:end))';
        if isfield(trajs,'U')
            U = shortdata(Ntrajs,3,Ubroad);
            xfakep = [longdata(xfakeshort(:,:,1:end-1))';...
                longdata(U(:,:,1:end-1))';...
                ones(1,(size(xfakeshort,3)-1)*Ntrajs,'like',Ybroad)];
            betaTRANS = (xfakef*xfakep')/(xfakep*xfakep');
            LDSparams.A = betaTRANS(:,1:Nx) +...
                randn(Nx)/50;
            LDSparams.B = betaTRANS(:,(Nx+1):(Nx+Nu)) +...
                randn(Nx,Nu)/50;
        else
            xfakep = [longdata(xfakeshort(:,:,1:end-1))';...
                ones(1,(size(xfakeshort,3)-1)*Ntrajs,'like',Ybroad)];
            betaTRANS = (xfakef*xfakep')/(xfakep*xfakep');
            LDSparams.A = betaTRANS(:,1:end-1);
        end
        LDSparams.muX = betaTRANS(:,end);
        LDSparams.SigmaX = (xfakef*xfakef') - betaTRANS*(xfakep*xfakef');
        
        % init ICs based on initial observations
        foo = randn(Nx,Nx,'like',Ybroad)/1000;         % small & random
        LDSparams.Info0 = pinv(cov((pinv(C)*Y0)')) + foo'*foo; % so it's invertible!
        LDSparams.mu0 = mean(pinv(C)*Y0,2);
        
        % tell ML4LDS that muX may be nonzero (see NOTE below)
        %%%
        LDSparams.muX = zeros(Nx,1,'like',Ybroad);
        LDSparams.muYX = zeros(Ny,1,'like',Ybroad);
        %%%
        
        % make the output and transition noises small [wrt what??]
        foo = randn(Ny,'like',Ybroad)/10;
        LDSparams.SigmaYX = foo'*foo;
        %%% NB: this will be overwritten if LDSdata has a field SigmaYX


    case 'random'
        
        LDSparams.A = randn(Nx,Nx,'like',Ybroad);
        if isfield(trajs,'U'), LDSparams.B = randn(Nx,Nu,'like',Ybroad); end
        foo = randn(Nx,Nx,'like',Ybroad);
        LDSparams.SigmaX = foo'*foo; % 1.0e-05*[0.5 0; 0 0.5];
        LDSparams.muX = randn(Nx,1,'like',Ybroad); % [0;0];
        LDSparams.C = randn(Ny,Nx,'like',Ybroad); % [1 0];
        foo = randn(Nx,Nx,'like',Ybroad);
        LDSparams.Info0 = foo'*foo; % inv([0.0030 0; 0 0.1000]);
        LDSparams.mu0 = randn(Nx,1,'like',Ybroad); % [1.0236; 2.5];
        if Ny > Nx, LDSparams.muYX = randn(Ny,1,'like',Ybroad); end
        
        % make the output and transition noises small [wrt what??]
        foo = randn(Ny,'like',Ybroad)/10;
        LDSparams.SigmaYX = foo'*foo;
        %%% NB: this will be overwritten if LDSdata has a field SigmaYX

   
    case 'FactorAnalysis'

        fprintf('Initializing LDS with factor analysis\n');
        [LDSparams.C,bhat,LDSparams.SigmaYX,What] = FactorAnalysis(Ybroad,Nx);
        Xhat = What*(Ybroad - bhat);
        
        Y = shortdata(Ntrajs,3,(Ybroad - bhat)');
        Z = shortdata(Ntrajs,3,Xhat');
        LDSparams = learnfullyobservedLDS(Y,Z);
        %%%%
        % But this will find a LDSparams.SigmaYX that is tiny!!!!
        %%%%

        if Ntrajs < length(LDSparams.mu0)
            LDSparams.mu0 = mean(Xhat,2);
            LDSparams.Info0 = inv(cov(Xhat')');
        end

end

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
[~,indices] = sort(diag(D),'descend');          % sort by eigenvalue
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
function [xpctT,XNtrpY] = Estep4LDS(trajs,LDSparams,figureInfo)

% booleans
BIASEDTRANSITIONS = 1;
BIASEDEMISSIONS = diff(size(LDSparams.C)) < 0; %%% see notes below

% Ns
Ntrajs = length(trajs);
Nstates = size(LDSparams.A,1);

% init
mucov2op = @(muA,muB,Sgm)(muA*muB' + sum(Sgm,3));
if isfield(trajs,'U')
    KFfunc = @(i,KFparams)(KalmanFilter(KFparams,trajs(i).Y,...
        'controls',trajs(i).U));
else
    KFfunc = @(i,KFparams)(KalmanFilter(KFparams,trajs(i).Y));
end
if ~isfield(trajs,'SigmaYX')
    [trajs(1:Ntrajs).SigmaYX] = deal(LDSparams.SigmaYX);
end
muX0 = 0; muXp = 0; muXf = 0; muX = 0; muYX = 0; XNtrpY = 0;
X0X0 = 0; XpXp = 0; XfXf = 0; XfXp = 0; XX = 0; YX = 0; YY = 0;
Nsamples = 0;

% loop through cases
for iTraj = 1:Ntrajs
    
    Nsamples = Nsamples + size(trajs(iTraj).Y,2);
    
    % filter
    LDSparams.SigmaYX = trajs(iTraj).SigmaYX;
    LDSparams.T = size(trajs(iTraj).Y,2);
    filtered = KFfunc(iTraj,LDSparams);
    
    
    % Approximate:
    %   Cov[X_{t+1},X_t| y_0,...,y_T] \approx Cov[X_{t+1},X_t| y_0,...,y_t] 
    %                                       = A*Cov[X_t|y_0,...,y_t]
    smoothed.XfXpCVRN = tensorOp(repmat(LDSparams.A,[1,1,LDSparams.T]),...
        filtered.CVRNMU);
    % WTF??  This appears totally wrong:
    %     % fake smoother
    %     for t = 1:(size(trajs(iTraj).Y,2)-1)
    %         Pt = filtered.CVRNMU(:,:,t);
    %         Jt = Pt*LDSparams.A'*filtered.INFOTU(:,:,t+1);
    %         XfXpCVRN = filtered.CVRNMU(:,:,t+1)*Jt'; %%% this is the approx.
    %         smoothed.XfXpCVRN(:,:,t) = Pt + Jt*XfXpCVRN - Jt*LDSparams.A*Pt;
    %         keyboard
    %     end
    
    
    smoothed.XHAT = filtered.XHATMU;
    smoothed.XCVRN = filtered.CVRNMU;
    
    
    % gather expected sufficient statistics
    % first-order
    muX0 = muX0 + smoothed.XHAT(:,1);
    muXp = muXp + sum(smoothed.XHAT(:,1:end-1),2);
    muXf = muXf + sum(smoothed.XHAT(:,2:end),2);
    muX  = muX  + sum(smoothed.XHAT,2);
    muYX  = muYX  + sum(trajs(iTraj).Y,2);
    
    % second-order
    X0X0 = X0X0 + mucov2op(smoothed.XHAT(:,1),smoothed.XHAT(:,1),smoothed.XCVRN(:,:,1));
    XpXp = XpXp + mucov2op(smoothed.XHAT(:,1:end-1),smoothed.XHAT(:,1:end-1),smoothed.XCVRN(:,:,1:end-1));
    XfXf = XfXf + mucov2op(smoothed.XHAT(:,2:end),smoothed.XHAT(:,2:end),smoothed.XCVRN(:,:,2:end));
    XfXp = XfXp + mucov2op(smoothed.XHAT(:,2:end),smoothed.XHAT(:,1:end-1),smoothed.XfXpCVRN);
    XX = XX + mucov2op(smoothed.XHAT,smoothed.XHAT,smoothed.XCVRN);
    YX = YX + mucov2op(trajs(iTraj).Y,smoothed.XHAT,0);
    YY = YY + mucov2op(trajs(iTraj).Y,trajs(iTraj).Y,0);
    
    % (average) log likelihood
    XNtrpY = XNtrpY + filtered.XNtrpY;
    
    if ~isempty(figureInfo)&&(iTraj==1)
        shatKF = LDSparams.C*filtered.XHATMU;
        shatRTSS = LDSparams.C*smoothed.XHAT;
        y = trajs(iTraj).Y;
        
        % now write to the plot
        %%%%% hardcoded fignum = 1
        animatePlot(figureInfo.num,figureInfo.subplotHandle,...
            figureInfo.dataplotHandle,1:size(y,2),y,shatKF,shatRTSS);
        % animatePlot(2,subplotHandleX,dataplotHandleX,1:T,zeros(ppp,T),... % x(1:ppp,:),...
        %     XhatKF(1:ppp,:),RTSSdstrbs.XHAT(1:ppp,:));
        % keyboard
        % pause();
        % end
    end
end


% now renormalize by the number of samples in each
XNtrpY = XNtrpY/Ntrajs;
xpctT.mu0 = muX0/Ntrajs;
xpctT.x0x0 = X0X0/Ntrajs;
xpctT.XpXp = XpXp/(Nsamples-Ntrajs);
xpctT.XfXf = XfXf/(Nsamples-Ntrajs);
xpctT.XfXp = XfXp/(Nsamples-Ntrajs);
if BIASEDTRANSITIONS
    xpctT.XpXp = [xpctT.XpXp muXp/(Nsamples-Ntrajs); muXp'/(Nsamples-Ntrajs) 1];
    xpctT.XfXp = [xpctT.XfXp muXf/(Nsamples-Ntrajs)];
end
xpctT.XX = XX/Nsamples;
xpctT.YX = YX/Nsamples;
xpctT.YY = YY/Nsamples;
if BIASEDEMISSIONS
    xpctT.XX = [xpctT.XX muX/Nsamples; muX'/Nsamples 1];
    xpctT.YX = [xpctT.YX muYX/Nsamples];
end



end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function animatePlot(fignum,subplotHandle,dataplotHandle,xdata,varargin)

% params
%%% Ndatavecs = length(varargin);
[Ndims,Ndatavecs] = size(dataplotHandle);

% pull up the current figure
set(0,'CurrentFigure',figure(fignum))
colors = 'krbcmy';
%%% try not to use any more than that

for ivec = 1:Ndatavecs
    
    ydata = varargin{ivec};
    %%%Ndims = size(ydata,1);
    
    for iDim = 1:Ndims
        set(dataplotHandle(iDim,ivec),'XData',xdata,'YData',ydata(iDim,:),...
            'color',colors(ivec));
    end
end
%%% set(subplotHandle,'XLim',[0 Max]);
%%% if you want to keep the x or y limits fixed over the animation....

end
%-------------------------------------------------------------------------%




%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Retired (for now, at least)
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
% There's a redundancy here.  A translation of Y (muYX!=0) can usually be
% absorbed into a constant input to X, namely muX := (I-A)C'inv(CC')muYX.
% Proceeding backwards from the conclusion:
%
%       Z_{t+1} = A*Z_t + (I-A)*C'*inv(C*C')*muYX
%       Y_t = C*Z_t
%   =>  Z_{t+1} = A*(Z_t - C'*inv(C*C')*muYX) + C'*inv(C*C')*muYX
%   =>  Z_{t+1} - C'*inv(C*C')*muYX = A*(Z_t - C'*inv(C*C')*muYX)
%   X := Z - C'*inv(C*C')*muYX =>
%       X_{t+1} = A*X_t.
%       Y_t = C*(X_t + C'*inv(C*C')*muYX)
%           = C*X_t + muYX.
%
% NB: *But this assumes that C*C' is invertible* ("C is fat").
%
% Likewise, a constant input to X (muX!=0) can be absorbed into a
% translation of the emission Y, namely muYX = C*inv(I-A)*muX:
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
BIASEDEMISSIONS = 0;    % isfield(LDSparams,'muYX');


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
