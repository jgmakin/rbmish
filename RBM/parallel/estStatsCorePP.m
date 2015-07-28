function varargout = estStatsCorePP(S,params,ESTIMATOR,varargin)
% ESTSTATSCOREPP    Core error stats code, parallel processing
%   ESTSTATSCOREPP turns data---population activities---in varargin into
%   errors by decoding and taking the difference with the underlying
%   stimuli in x via the subroutine estError.m.  It needs also needs to
%   know what ESTIMATOR to use, and the params data structure.
%
%  NB: D and S are expected in longdata rather shortdata form; see 
%  longdata.m

%-------------------------------------------------------------------------%
% Revised: 06/16/14
%   -Nmods = params.Nmods -> Nmods = size(S,3).  This should almost never
%   make any difference, but it seems more logical.
% Revised: 12/11/13
%   -changed input argument x to S, which has an extra dimension (two dims
%   of params.Ndims and params.Nmods, rather than one of Ndims*Nmods); made 
%   the appropriate changes in the function to account for this.
% Revised: 12/10/13
%   -changed some indexing (but the results are still the same)
%   -changed the form of all the outputs!!  E.g., the stats structures were
%   formerly in a cell array: s1in, s2in, s3in, s1out, s2out, s3out.  Now
%   they're in a cell matrix, {s1in, s1out; s2in s2out; ...}
%   -shatL and shatN are also differently shaped, b/c estError.m has
%   changed.
% Revised: 05/05/11
%   -changed it to varargin so that it might accomodate only one dataset
%   rather than necessarily both input and output.
% Adapted: 05/04/11
%   -from estimatorStatsPP
%-------------------------------------------------------------------------%
% Revised: 01/31/11
%   -changed multidimensional array e from 4 to 3 dims
%   -changed to single ouput (since all errors are now computed in their
%   local modalities; see estError.m)
% Revised: 12/20/10
%   -changed to accomodate 3-modality scheme
% Revised: 11/08/10
%   -reshaped D to precompute the outputs in one go instead of in a
%   loop
% Revised: 09/08/10
%   -calls the new estError.m [see] twice rather than once
% Created: 07/07/10
%   by JGM
%-------------------------------------------------------------------------%

% init
n = nargin - 3;
Nmods = size(S,3);
[pool,HADBEENCLOSED] = parallelInit;

% malloc
statsL = cell(Nmods,n);
statsN = cell(Nmods,n);

% loop
for i = 1:n
    D = varargin{i};
    Nexamples = size(D,1);
    parfor (j = 1:Nexamples, 8)
        s = shiftdim(S(j,:,:),1);
        [eL(j,:,:),eN(j,:,:),shatL(j,:,:,i),shatN(j,:,:,i)] =...
            estError(D(j,:),s,params,ESTIMATOR,0);
    end
    for iMod = 1:Nmods
        statsL{iMod,i}.cov = nancov(eL(:,:,iMod));    
        statsL{iMod,i}.mu = nanmean(eL(:,:,iMod))';
        statsN{iMod,i}.cov = nancov(eN(:,:,iMod));    
        statsN{iMod,i}.mu = nanmean(eN(:,:,iMod))';
        %%% maybe it's a little strange to have your data stored with
        %%% examples in the rows, and then the means stored as column
        %%% vectors....
    end
end

% close the pool, if it started closed
if HADBEENCLOSED, delete(pool); end

% *** (2) ***

% collect outputs
varargout{1} = statsL;
varargout{2} = statsN;
varargout{3} = shatL;
varargout{4} = shatN;


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% *** (1) ***
% For some insane reason, nargout changes its value within the parfor loop.
% I'd call this a bug.  So you have to assign it to some constant at init.
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
% *** (2) ***
% find the indices for which the Jacobian is a good approx., and keep those
% if nargout == 5
%     indexvec = errorPrune(eLocal,eNeutr,x,params);
%     eLocal = eLocal(:,:,indexvec);
%     eNeutr = eNeutr(:,:,indexvec);
%     shatL = shatL(indexvec,:);
%     shatN = shatN(indexvec,:);
% end
% ...
% if nargout == 5
%     varargout{5} = x(indexvec,:);
% end
%-------------------------------------------------------------------------%