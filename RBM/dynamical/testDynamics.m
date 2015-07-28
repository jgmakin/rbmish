function eStats = testDynamics(LDSdata,params,TOPLOT,varargin)
% testDynamics      Compute/plot error statistics for.....
%
% USAGE:
%   eStats = testDynamics(LDSdata,params,1,pSENSORY)
%   eStats = testDynamics(LDSdata,params,1,pSENSORY,pEFH,pKFobs,pKFem)

%-------------------------------------------------------------------------%
% Revised: 03/30/15
%   -fixed bug for assigning modality names to eStats.tags.name
% Revised: 09/14/14
%   -added eStats as an output argument
% Created: 07/15/14
%   by JGM
%-------------------------------------------------------------------------%


% use unwrapped stimuli for case 'wrapping'
if strcmp(params.dynamics.walls,'wrapping'),
    S = unwrapStims(LDSdata,params);
else
    S = LDSdata.S;
end
S = longdata(S);

% get the expectations and names
[Nexamples,Ndims,Nmods] = size(S);
Nposteriors = length(varargin);
posteriorArray = [varargin{:}];
allXpct = longdata(cat(2,posteriorArray(:).Xpct));
allXpct = reshape(allXpct,[Nexamples,Ndims,Nmods,Nposteriors]);
allNames = {posteriorArray.name};


% loop through modalities (probably PROP and CTRL)
for iMod = 1:Nmods
    data.Xpct(1:Nexamples,1:Ndims,1:Nposteriors) = allXpct(:,:,iMod,:);
    data.srcs = allNames;
    if sum(strcmp('sensory',allNames))
        [data.srcs{strcmp('sensory',allNames)}] = deal(params.mods{iMod});
    end
    eStats(iMod) = getErrorStats(data,S(:,:,iMod));
    if TOPLOT
        fprintf('\nwriting tikz plot for %s...\n',params.mods{iMod});
        dispErrStats(eStats(iMod),params.mods{iMod});
    else
        fprintf('\nskipping tikz plots...\n');
    end
end



end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function S = unwrapStims(ldsDATA,params)
% it's more like "use unwrapped versions" than "unwrap"

if isfield(params.dynamics,'H') %%% not ideal
    C = blkdiag(params.dynamics.C,params.dynamics.H);
else
    C = params.dynamics.C;
end
S = reshape(shortdata(params.Ncases,3,longdata(ldsDATA.Z)*C'),...
    size(ldsDATA.S));


end
%-------------------------------------------------------------------------%
