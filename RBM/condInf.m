function [statsOnePass,statsFinal] = condInf(D0,S0,wts,deadstim,params)
% compute error statistics for "conditional inference," i.e. inferring
% visual hand position from proprioceptive hand position, and vice versa.
% This operation is sometimes referred to as "function approximation."
%
% USAGE:    load wtsfile.m
%           [D0,S0] = generateData(500,params);
%           [statsOnePass statsFinal] = condInf(D0,x0,wts,deadstim,params)
% 
% NB: THIS FUNCTION IS CURRENTLY BROKEN!

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -X -> S and corresponding changes
% Revised: 12/12/13
%   -changed to accommodate new orientation of Dlong
% Revised: 08/16/12
%   -rewrote to accomodate length(params.mods)=3.
%   -replaced while loop with recursive updateInputs.
%   -etc.
% Revised: 08/17/11
%   -incorporated "sampling with replacement" (constant s, varying r)
% Adapted: 8/16/10
%   -from the getNoiseFloor.m and unimodalstim.m
%   by JGM
%-------------------------------------------------------------------------%


%%%%%%%%%%%%%%
%%%% This function was broken by changes to other code, e.g. the new
%%%% version of PPCencode.  You should re-write it.
%%%%%%%%%%%%%%


% init
params.resample = 1;

% you have to process the examples severally, so no use having batches
[Di,Si] = longdata(D0,S0);
clear D0 S0;
[Nexamples,Ndims] = size(Di);

% zero out one input modality, for all examples
[Di,S.zvec] = killoneinput(Di,deadstim,length(params.mods));


% stats for a single up-down pass (use neutral-space stats only)
Do = updownDBN(Di,wts,params,'suffstats');
[~, statsOnePass] = estStatsCorePP(Si,params,'CoM',Di,Do);

% compute one-pass (normalized) reconstruction error
E = Di - Do;
reconErr = sqrt(sum(E.*E,2))/Ndims;


% malloc (unnecessary if using parfor)
Dfinal = zeros(size(Di));

% matlabpool open local 8
% parfor (iExample = 1:Nexamples, 8)
for iExample = 1:Nexamples
    
    [S.xPatch,S.gains] = getWorldVars(squeeze(Si(iExample,:,:)),params);
    dClamped = updateWithClamp(Di(iExample,:),Do(iExample,:),S,params); 
    [Dfinal(iExample,:),counter] =...
        updateInputs(dClamped,reconErr(iExample,:),1,S,wts,params);
    fprintf('counter exit: %i  example: %i\n',counter,iExample);
    
end
% matlabpool close

% stats for the "converged" output (use neutral-space stats only)
[~, statsFinal] = estStatsCorePP(Si,params,'CoM',Di,Dfinal);
%%%%%% seems like you could leave out Di from this function...


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [Di,zvec] = killoneinput(Di,deadstim,Nmods)

% init
Ndims = size(Di,2);

% set one modality to zero for all examples
switch deadstim
    case 'vis',     zvec = 1:Ndims/Nmods;
    case 'prop',    zvec = (Ndims/Nmods+1):(2*Ndims/Nmods);
    case 'eye',     zvec = (Ndims-Ndims/Nmods+1):Ndims;
    otherwise, error('unrecognized stimulus (should be ''vis,'' ''prop,'' or ''eye''');
end
Di(:,zvec) = 0;


end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [d_out,counter] = updateInputs(d_in,eOLD,counter,S,wts,params)

% init
alpha = 0.001; % /100000/100;

% increment the counter and compute the new error
counter = counter + 1;
d_out = updownDBN(d_in,wts,params,'means','quiet');
% [~, d_out] = updownfast(d_in,wts,params);
eNEW = norm(d_in - d_out)/length(d_in);

% stop iterating?
if (eOLD - eNEW)/eOLD < alpha % &&0
    % if change in error is small
    % fprintf('error (eNew = %f, eOLD = %f) stopped decreasing\n',eNEW,eOLD);
    return;
elseif counter > 100
    % don't wait forever for convergence
    fprintf('counter ran out!\n');
    return;
else
    
    % draw a new set of spike counts for the same S
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
Nmods = length(params.mods);
Nvis = length(dClamped);

% resample spike counts (stim fixed)?
if params.resample
    
    d = zeros(Nvis/Nmods,Nmods);
    for j = 1:Nmods
        d(:,j) = PPCencode(S.xPatch(:,j),S.gains(j),...
            params.typeUnits{1},params.numsUnits{1},params);
    end
    dClamped = d(:);
end

% but update the modality that you're inferring
dClamped(S.zvec) = d_out(S.zvec);


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [xPatch,gains] = getWorldVars(s,params)

% init
Nmods = length(params.mods);
Ndims = params.Ndims;

% compute useful quantities for replacement sampling
if params.resample
    smin = params.smin;
    smax = params.smax;
    patchmin = params.margin*ones(1,Ndims);
    patchmax = patchmin + params.respLength;   
    gains = mean([params.gmin; params.gmax]);
    
    xPatch = zeros(Ndims,Nmods);
    for j = 1:Nmods
        xPatch(:,j) = scalefxn(s(:,j),smin(:,j),smax(:,j),patchmin,patchmax);
    end
end
%%%%
% need an else clause for this to make any sense---otherwise you could have
% unassigned outputs.  If you *know* that params.resample=1, don't bother
% to put in the original if clause.
%%%%

%%%%
% should this new gain really be fixed, or should it be samples?
%%%%

end
%-------------------------------------------------------------------------%

