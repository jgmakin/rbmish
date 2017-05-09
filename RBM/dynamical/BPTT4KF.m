function LDSparams = BPTT4KF(YSTAR,LDSdata,LDSparams,mAC)
% BPTT4KF   "Backpropagation through time" for the Kalman-filter equations
%
% USAGE:
%   LDSparams = BPTT4KF(YSTARtrain,pMirrorTrain.XHATMU,...
%        pMirrorTrain.CVRNMU,LDSdataTrain,LDSparams,2e-5/T);
%
% See fitLDStoKFobsvEsts.m, which calls this function.  This function
% implements gradient descent on a loss function defined in terms of the
% Kalman-filter recursions.  In particular, it penalizes squared deviation
% of a KF's predicted emission, yhat = C*zhat, from a desired predicted
% emission, ystar.  Doing so requires a backwards propagation of
% derivatives; these are computed in your lab notes.
%
% The last argument is the inverse learning rate.


%-------------------------------------------------------------------------%
% Revised: 06/25/15
%   -made parameter updates happen every *trajectory*, rather than every
%   data set; now converges much faster.
% Created: 06/24/15
%   by JGM
%-------------------------------------------------------------------------%


% extract the parameters
A = LDSparams.A;
C = LDSparams.C;
SigZ = LDSparams.SigmaX;
% let's just assume that muX = 0...

% Ns
[Ntraj,~,T] = size(YSTAR);
[Ny,Nz] = size(C);

% params
TOPLOT = 0;

% malloc
g = zeros(Nz,T,'like',YSTAR);
B = zeros(Nz,Nz,T,'like',YSTAR);
F = zeros(Nz,Nz,T,'like',YSTAR);
Qinv = zeros(Ny,Ny,T,'like',YSTAR);
Wtut = zeros(Nz,Nz,T,'like',YSTAR);
Wyt = zeros(Nz,Ny,T,'like',YSTAR);


% loop through trajectories
fprintf('Learning parameters')
for iTraj = 1:Ntraj
    
    if any(eig(SigZ)<0)
        eig(SigZ)
        error('uh-oh: SigmaZ is not positive definite!');
    end
    fprintf('.');

    
    % Kalman filter this trajectory
    y = shiftdim(LDSdata.Y(iTraj,:,:),1);
    SigYZ = shiftdim(LDSdata.SigmaYX(iTraj,:,:,:),1);
    KF = KalmanFilter(setfield(LDSparams,'SigmaYX',SigYZ),y);
    zhat = KF.XHATMU;
    P = KF.CVRNMU;
    
    % collect errors for this trajectory
    err = C*zhat - shiftdim(YSTAR(iTraj,:,:),1);
    %%% fprintf('traj: %i; trerror: %2.4d\n',iTraj,mean(err.^2));
    
    % loop backwards through time
    for t = (T-1):-1:1
        
        % update *future* params
        Rtplus1 = A*P(:,:,t)*A' + SigZ;
        Qinv(:,:,t+1) = inv(SigYZ(:,:,t+1) + C*Rtplus1*C');
        Wyt(:,:,t+1) = Rtplus1*C'*Qinv(:,:,t+1);
        Wtut(:,:,t+1) = P(:,:,t+1)/Rtplus1;
        F(:,:,t+1) = Wtut(:,:,t+1)'*(g(:,t+1)*(zhat(:,t+1) -...
            A*zhat(:,t))'/P(:,:,t+1) + B(:,:,t+1))*Wtut(:,:,t+1);
        
        % the backwards recursions
        g(:,t) = 2*C'*err(:,t) + A'*Wtut(:,:,t+1)'*g(:,t+1);
        B(:,:,t) = A'*F(:,:,t+1)*A;
        
    end
    
    if TOPLOT
        figure(45); clf; hold on;
        plot(g');
        hold off;
        drawnow;
    end
    
    
    % assemble pieces for parameter updates
    FplusF = F + permute(F,[2,1,3]);                        % Nz x Nz x T
    FA = tensorOp(FplusF,repmat(A,[1,1,T]));                % Nz x Nz x T
    FAP = tensorOp(FA(:,:,2:end),P(:,:,1:end-1),[1,3,2]);   % Nz x Nz x T-1
    Wtutg = permute(tensorOp(permute(Wtut,[2,1,3]),...
        permute(g,[1,3,2])),[1,3,2]);                       % Nz x T
    
    delA = (Wtutg(:,2:end)*zhat(:,1:end-1)' + sum(FAP,3))/T;% Nz x Nz
    delSigZ = (sum(FplusF,3) - diag(diag(sum(F,3))))/T;     % Nz x Nz
    
    % pretty unfortunate
    diffQinv = permute(tensorOp(...
        permute(C*A*zhat(:,1:end-1) - y(:,2:end),[3,1,2]),...
        Qinv(:,:,2:end)),[2,3,1]);                          % Ny x T-1
    Wytgt = permute(tensorOp(permute(Wyt,[2,1,3]),...
        permute(g,[1,3,2])),[1,3,2]);                       % Ny x T
    WytBt = tensorOp(permute(Wyt,[2,1,3]),B);               % Ny x Nz x T
    
    delSigYZ = (Wytgt(:,2:end)*diffQinv' +...
        sum(tensorOp(WytBt(:,:,2:end),Wyt(:,:,2:end)),3))/T;% Ny x Ny
    
    gtPt = permute(tensorOp(permute(g,[3,1,2]),P),[2,3,1]); % Nz x T
    BBP = tensorOp(B + permute(B,[2,1,3]),P);               % Nz x Nz x T
    foo1 = 2*err(:,2:end)*zhat(:,2:end)';                   % Ny x Nz
    foo2 = -diffQinv*gtPt(:,2:end)';                        % Ny x Nz
    foo3 = -Wytgt(:,2:end)*zhat(:,2:end)';                  % Ny x Nz
    foo4 = -sum(tensorOp(permute(Wyt(:,:,2:end),[2,1,3]),...
        BBP(:,:,2:end)),3);                                 % Ny x Nz
    delC = (foo1 + foo2 + foo3 + foo4)/T;                   % Ny x Nz
    
    
    % update params
    A = A - delA/mAC;
    C = C - delC/mAC;
    SigZ = SigZ - delSigZ/mAC^2;
    %%% It appears necessary either to delay the modification of SigZ till 
    %%% after a decent A and C have been learned, or to make its learning
    %%% rate much smaller.  The use of a square seems intuitive, since SigZ
    %%% is in squared units (compare the standard deviations).
    % delSigYZ = delSigYZ + delSigYZ' - diag(diag(delSigYZ));
    % SigYZ = SigYZ - delSigYZ/mw;
    LDSparams.A = A;
    LDSparams.C = C;
    LDSparams.SigmaX = SigZ;
end
fprintf('\n');










