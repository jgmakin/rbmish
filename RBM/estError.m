function varargout = estError(d,strueL,params,ESTIMATOR,DISP)
% ESTERROR  MLE error for PPC models
%   USAGE: [eL eN] = estError(d,s,params,'CoM',0)
%          [eL eN shatL shatN] = estError(d,x,params,'CoM',0)
%
%   ESTERROR computes the error of estimators for an underlying stimulus
%   (strueL).  The estimators are derived from population activities d
%   according to the method prescribed in the string ESTIMATOR.
%
%   d and strueL are structured in "the usual way" (see the rest of yr
%   RBM code, esp. DATAGEN.m and estimatorStats.m).  d is assumed to
%   contain a *pair* of populations, one for each transducer (vis, prop),
%   and likewise for strue (*in its columns*)
%
%   For proper comparison, the prop estimates are transformed into vis
%   estimates (and vice versa) and stored in a "visual error" (e_v)
%   variable (and likewise for e_p).  The variable e is a "mixed error"
%   (i.e., does not apply the xformations), and is probably useless.
%
%   DISP optionally toggles display of the population activities, overlayed
%   with the true underlying stimulus as well as the transformed (prop to
%   vis) estimate.
%
%   The output e is a multidimensional array:
%
%       row:    dimension of the stimulus; max = params.m
%       column: input (vis, prop); max = 2
%       slice:  1 (used for iCase in e.g. estimatorStats.m)
%       block:  which space/units (vis, prop, eye, etc); max = params.Nmods


%-------------------------------------------------------------------------%
% Revised: 12/10/13
%   -changed input from xtrue to strue (different dimensions)
% Revised: 12/10/13
%   -changed outputs shatL and shatN to be matrices (Ndims x Nmods) rather than
%   vectors.
% Revised: 04/08/11
%   -reverted yet again to to computing the errors in both local and
%   foreign ("neutral") spaces; cf. labnotes.pdf
% Revised: 01/31/11
%   -eliminated estGather; now all errors are computed in their local
%   modality and not in a common mod.
% Revised: 01/24/11
%   -added selfmark.m and varargin{2} for display mode
% Revised: 12/20/10
%   -changed to accomodate 3-modality scheme
% Revised: 11/29/10
%   -changed to accomodate x in *true* coordinates
%   -replaced <estimate> with decode.m
% Revised: 09/08/10
%   -rewrote to operate on a single input dataset, D (so now estimatorStats
%   has to call it twice)
%   -eliminated "redundancy" in estimate.m (and etc.) by making it run on a
%   single half-population alone
% Cribbed: 08/27/10
%   from estimatorStats
%   by JGM
%-------------------------------------------------------------------------%

% init
Nmods = params.Nmods;
Ndims = params.Ndims;
smin = params.smin;
smax = params.smax;
T = displayshape(d,params);                     % vec -> topo.


% compute estimators from PPCs
% local
shatL = zeros(Ndims,Nmods);
for iMod = 1:Nmods
    shatL(:,iMod) = decode(T{iMod},[smin(:,iMod) smax(:,iMod)],params,ESTIMATOR);
end

% "neutral"
strueN = estGather(strueL,params);
shatN = estGather(shatL,params);

% compute error
eL = shatL - strueL;
eN = shatN - strueN;


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this doesn't look very robust (jgm, 12/10/13)
% plot?
if DISP
    PPCplot(cat(2,T{:}),params,'(vis, prop)');
    hold on;
    selfmark(shatL,strueL,max(cat(2,T{:})),params);
    hold off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%

varargout{1} = eL;
varargout{2} = eN;
varargout{3} = shatL;
varargout{4} = shatN;


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function selfmark(shat,strue,rmax,params)
% SELFMARK  Mark populations with the estimates they encode
%   SELFMARK takes a (displayshaped) representation T of a population code
%   for some modality, along with the true stimulus xtrue that the pop.
%   encodes, and plots the latter along with the the decoded estimates on
%   top of the former.

%%%%%%%%%%%%%%%%
% Broken, at least b/c of new shape of s (jgm, 12/16/13)
%%%%%%%%%%%%%%%%


% params
Nmods = params.Nmods;
smin = params.smin;
smax = params.smax;

% loop through vis and prop
s = cat(3,shat,strue);
c = 'rg';
for j = 1:Nmods
    for i = 1:2                     % first the ests, then the true stimuli
        mark(s(:,j,i),[smin(:,j),smax(:,j)],rmax(j),j-1,c(i),params);
    end
end



end
%-------------------------------------------------------------------------%