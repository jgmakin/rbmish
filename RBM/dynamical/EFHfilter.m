function varargout = EFHfilter(LDSdata,wts,params)
% USAGE:
%
%   [V0, pSENSORY] = EFHfilter(LDSdata,wts,params)
%   [V0, pSENSORY, pRBM] = EFHfilter(LDSdata,wts,params)
%   [V0, pSENSORY, pRBM1, pRBM2] = EFHfilter(LDSdata,wts,params)
%
% Each pRBM iterates through one further step of Gibbs "meaning."
%
% What is stored in p?  Well,
%
%       p(s|r) = N(CoM,tuningCov/eta)       (assuming flat prior)
%   =>  p(x|r) = N(pinv(C)*CoM,"pinv(C)*tuningCov/eta*pinv(C)")'
%   =>  p(C*x|r) = N(CoM,tuningCov/eta)
%
% The posterior mean E[X|r] =: Xhat; the posterior mean E[C*X|r] =: CXhat.
% The quantity "pinv(C)*tuningCov/eta*pinv(C)" isn't really the posterior
% covariance Cov[X|r] b/c of the pseudo-inverse...  One could instead store
% the information matrix....  However we do store Cov[C*x|r] =: CvrnCXhat.


%-------------------------------------------------------------------------%
% Revised: 01/29/15
%   -Ncases is now set by the size of LDSdata, rather than params.Ncases
% Revised: 01/03/14
%   -made it actually work with multiple "modalities" (probably prop and
%   the control)
% Revised: 11/05/13
%   -changed the input and output arguments!!, rationalized
% Revised: 07/03/13
%   -changed the way the RESTART flag works
%   -changed (back) to varargout, so as to pass out the filterdata
% Revised: 07/01/13
%   -added force-to-zero for the restarts (in case of dynamics.walls =
%   'restart')
%   -added varargout for the hidden layer activities
% Cribbed: 07/01/13
%   from getFilteredData
%   by JGM
%-------------------------------------------------------------------------%

% N's
Nvis = params.numsUnits(1);
Ncases = size(LDSdata.R,1);
Niter = nargout-2;
T = size(LDSdata.R,3);

% other params
vtype = params.typeUnits{1};
htype = params.typeUnits{2};
switch params.typeUnits{1}
    case {'BP','GP','Bernoulli'}, inds = (params.t+1):Nvis;
    case 'PB', inds = 1:(Nvis-params.t);
    otherwise, error('strange type of units for this function -- jgm');
end


% malloc
R = zeros([size(LDSdata.R),Niter],'like',LDSdata.R);
thisV0 = zeros(Ncases,params.numsUnits(2),'like',LDSdata.R);
V0 = zeros(Ncases,params.numsUnits(2),T,'like',LDSdata.R);
%%% thisR1 = zeros(Ncases,Nvis,'like',wts);

% up-down pass
fprintf('\nFiltering with EFH');
for t = 1:T
    
    % get the visible (0), hidden, and visible (1) activities
    thisR0 = LDSdata.R(:,:,t);
    
    if Niter > 0
        thisRecurrence = thisV0;
        %%% thisRecurrence = thisR1(:,1:Nhid); % assumes left-half...
        
        % force recurrent activities to 0 upon restarts
        thisRecurrence(LDSdata.restarts{t},:) = 0;
        
        % all visible units
        thisY = [thisRecurrence, thisR0];
        thisV0 = feedforward(thisY,wts{1}(1:end-1,:),wts{1}(end,:),htype,params);
        
        % store
        V0(:,:,t) = thisV0;
    end
    
    % Gibbs ``meaning''
    for iIter = 1:Niter
        
        % one Gibbs step
        thisV = feedforward(thisY,wts{1}(1:end-1,:),wts{1}(end,:),htype,params);
        thisY = feedforward(thisV,wts{2}(1:end-1,:),wts{2}(end,:),vtype,params);
        
        % store the RHS
        R(:,:,t,iIter) = thisY(:,inds);
    end
    
    fprintf('.');
end
fprintf('\n');



% clear some space
clear thisRecurrence thisV thisY thisR thesePPCs D0

% sufficient stats from the naive decoder
p{1}.Xpct = LDSdata.Y;
p{1}.Cvrn = LDSdata.SigmaY;
p{1}.name = 'sensory';
p{1}.params = []; %%% maybe one day you'll use this

% sufficient stats from the rEFH (perhaps multiple iterations)
for iIter = 1:Niter
    [p{iIter+1}.Xpct,p{iIter+1}.Cvrn] =...
        PPC2FltrCmlnts(R(:,:,:,iIter),LDSdata.S,LDSdata.Z,params);
    p{iIter+1}.name = 'rEFH';
    p{iIter+1}.params = []; %%% maybe one day you'll use this...
end

% store these posteriors
varargout = {V0, p{:}};


end














