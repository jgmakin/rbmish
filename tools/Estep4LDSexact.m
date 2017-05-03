function [xpctT,XNtrpY,infrs] = Estep4LDSexact(trajs,LDSparams,figureInfo)

% booleans
BIASEDTRANSITIONS = 1;
BIASEDEMISSIONS = diff(size(LDSparams.C)) < 0; %%% see notes below

% Ns
Ntrajs = length(trajs);

% init
mucov2op = @(muA,muB,Sgm)(muA*muB' + sum(Sgm,3));
if isfield(trajs,'U')
    KFfunc = @(i,KFparams)(KalmanFilter(KFparams,trajs(i).Y,...
        'controls',trajs(i).U));
else
    KFfunc = @(i,KFparams)(KalmanFilter(KFparams,trajs(i).Y));
end
if ~isfield(trajs,'SigmaY')
    [trajs(1:Ntrajs).SigmaY] = deal(LDSparams.SigmaY);
end
if nargout > 2, infrs(Ntrajs,1).Z = []; end
muX0 = 0; muXp = 0; muXf = 0; muX = 0; muY = 0; XNtrpY = 0;
X0X0 = 0; XpXp = 0; XfXf = 0; XfXp = 0; XX = 0; YX = 0; YY = 0;
Nsamples = 0;

% loop through cases
for iTraj = 1:Ntrajs
    
    Nsamples = Nsamples + size(trajs(iTraj).Y,2);
    
    % filter and smooth
    LDSparams.SigmaY = trajs(iTraj).SigmaY;
    LDSparams.T = size(trajs(iTraj).Y,2);
    filtered = KFfunc(iTraj,LDSparams);
    smoothed = RTSsmoother(LDSparams,filtered);
    %%%foo(iCase).X = RTSS.XHAT;
    
    % gather expected sufficient statistics
    % first-order
    muX0 = muX0 + smoothed.XHAT(:,1);
    muXp = muXp + sum(smoothed.XHAT(:,1:end-1),2);
    muXf = muXf + sum(smoothed.XHAT(:,2:end),2);
    muX  = muX  + sum(smoothed.XHAT,2);
    muY  = muY  + sum(trajs(iTraj).Y,2);
    
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
    
    if nargout > 2, infrs(iTraj).Z = smoothed.XHAT; end
        
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
    xpctT.YX = [xpctT.YX muY/Nsamples];
end



end