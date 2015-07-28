function [statsOnePass,statsFinal] = condInfPP(D0,S0,wts,deadstim,params)
% compute error statistics for "conditional inference," i.e. inferring
% visual hand position from proprioceptive hand position, and vice versa.
% This operation is sometimes referred to as "function approximation."
%
% USAGE:    load wtsfile.m
%           [D0,S0] = DATAGENPP(500,params);
%           [statsOnePass statsFinal] = condInfPP(D0,S0,wts,deadstim,params)

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -X -> S and corresponding changes
% Revised: 12/12/13
%   -TO DO: dclamped needs to become a *column* vector....
% Revised: 08/16/12
%   -rewrote to accomodate params.Nmods=3.
%   -replaced while loop with recursive updateInputs.
%   -etc.
% Revised: 08/17/11
%   -incorporated "sampling with replacement" (constant s, varying r)
% Adapted: 8/16/10
%   -from the getNoiseFloor.m and unimodalstim.m
%   by JGM
%-------------------------------------------------------------------------%

% init
params.resample = 0;
VERBOSE = 1;

% you have to process the examples severally, so no use having batches
[Di,Si] = longdata(D0,S0);
clear D0 S0;
[Nexamples,Ndims] = size(Di);

% zero out one input modality, for all examples
[Di,S.zvec] = killoneinput(Di,deadstim,params.Nmods);


% stats for a single up-down pass (use neutral-space stats only)
[~, Do] = updown(Di,wts,params,'means','quiet');
[~, statsOnePass] = estStatsCorePP(Si,params,'CoM',Di,Do);

% compute one-pass (normalized) reconstruction error
E = Di - Do;
reconErr = sqrt(sum(E.*E,2))/Ndims;


% malloc (unnecessary if using parfor)
Dfinal = zeros(size(Di));

% open parallel pool
[pool,HADBEENCLOSED] = parallelInit;

parfor (iExample = 1:Nexamples, 8)
    T = getWorldVars(squeeze(Si(iExample,:,:)),S,params);
    dClamped = updateWithClamp(Di(iExample,:),Do(iExample,:),T,params);
    [Dfinal(iExample,:),counter] =...
        updateInputs(dClamped,reconErr(iExample),1,T,wts,params);
    if VERBOSE
        fprintf('counter exit: %i  example: %i\n',counter,iExample);
    end
    
end
if HADBEENCLOSED, delete(pool), end

% stats for the "converged" output (use neutral-space stats only)
[~, statsFinal] = estStatsCorePP(Si,params,'CoM',Di,Dfinal);
%%%%%% Seems like you could leave Di out of this function!!


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [d_out,counter] = updateInputs(d_in,eOLD,counter,S,wts,params)

% init
alpha = 0.001; % /100000/100;

% increment the counter and compute the new error
counter = counter + 1;
[~, d_out] = updown(d_in,wts,params,'means','quiet');
% [~, d_out] = updownfast(d_in,wts,params);
%%%%%%%
eNEW = norm(d_in - d_out)/length(d_in);
%%% Isn't this already in updown??
%%%%%%%


% stop iterating?
if (eOLD - eNEW)/eOLD < alpha % &&0
    % if change in error is small
    % fprintf('error (eNew = %f, eOLD = %f) stopped decreasing\n',eNEW,eOLD);
    d_out = d_in;
    return;
elseif counter > 100
    % don't wait forever for convergence
    fprintf('counter ran out!\n');
    return;
else
    
    % draw a new set of spike counts for the same x
    dClamped = updateWithClamp(d_in,d_out,S,params);
    
    % recurse
    [d_out,counter] = updateInputs(dClamped,eNEW,counter,S,wts,params);
    
end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function dClamped = updateWithClamp(dClamped,d_out,S,params)
%%%%%%%%%%%
% This function was probably broken by changing GTrespfxn to work on
% matrices rather than vectors---need to check!
%%%%%%%%%%%

% init
Nmods = params.Nmods;
Nvis = length(dClamped);

% resample spike counts (stim fixed)?
if params.resample
    d = zeros(Nvis/Nmods,Nmods);
    for j = 1:Nmods
        d(:,j) =...
            PPCencode(S.xPatch(:,j),S.gains(j),params.typeUnits{1},params);
    end
    dClamped = d(:)';
end

% but update the modality that you're inferring
dClamped(S.zvec) = d_out(S.zvec);


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function S = getWorldVars(Si,S,params)

% init
Nmods = params.Nmods;
Ndims = params.Ndims;

% compute useful quantities for replacement sampling
if params.resample
    smin = params.smin;
    smax = params.smax;
    patchmin = params.margin*ones(1,Ndims);
    patchmax = patchmin + params.respLength;
    
    gains = params.g*ones(Nmods,1);
    if isfield(params,'gains'), gains = params.gains; end
    
    xPatch = zeros(Ndims,Nmods);
    for j = 1:Nmods
        xPatch(:,j) = scalefxn(Si(:,j),smin(:,j),smax(:,j),patchmin,patchmax);
    end
    
    S.xPatch = xPatch;
    S.gains = gains;
end
%%%%
% should this new gain really be fixed, or should it be samples?
%%%%

end
%-------------------------------------------------------------------------%