function [stVWDR,stWDR] = adaptByRule(D,shatL0,etaVWDR,etaWDR,N,params)
% This function assumes that the covariances, covV and covP, are in the
% same space!  Likewise, shft must be in that space as well.

% init
tuningCov = computetuningcovs(params);
Ndims = params.Ndims;
Nmods = length(params.mods);

% malloc
bias = zeros(Ndims,Nmods);
stVWDR = reshape(zeros(Nmods*Ndims,N),Ndims,Nmods,N);
stWDR = reshape(zeros(Nmods*Ndims,N),Ndims,Nmods,N);


% loop
sLvwdr = shatL0;                sNvwdr = estGather(sLvwdr,params);
sLwdr = shatL0;                 sNwdr = estGather(sLwdr,params);
%%%% change to while loop?
for i = 1:N
    
    % get next predicted step (based on current location)
    ind = ceil(size(D,1)*rand);
    [stepVWDR, ~, ~] = getPredictedStep(D(ind,:),tuningCov,bias,...
        sLvwdr,sNvwdr,params);
    [~, stepWDR, ~] = getPredictedStep(D(ind,:),tuningCov,bias,...
        sLwdr,sNwdr,params);
    
    % update estimates
    sNvwdr = sNvwdr + etaVWDR*stepVWDR;
    sNwdr = sNwdr + etaWDR*stepWDR;
    
    % convert back to local space for SICE
    sLvwdr = sNvwdr;    sLvwdr(:,1) = FK2link(sNvwdr(:,1),params.roboparams,1)';
    sLwdr = sNwdr;      sLwdr(:,1) = FK2link(sNwdr(:,1),params.roboparams,1)';
    
    % collect
    stVWDR(:,:,i) = etaVWDR*stepVWDR;
    stWDR(:,:,i) = etaWDR*stepWDR;
    
    if ~mod(i,1000)
        fprintf('i = %i\n',i);
    end
end

end
%-------------------------------------------------------------------------%
