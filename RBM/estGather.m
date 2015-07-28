function sN = estGather(s,params)
% ESTGATHER    Gather estimators for each modality
%   ESTGATHER assembles all estimators of the stimulus, in the space of each
%   modality, into the multidimensional array S:
%
%       row:    dimension of the stimulus; max = params.m
%       column: input (vis, prop); max = params.Nmods


%-------------------------------------------------------------------------%
% Revised: 12/17/13
%   -added case for MODEL 'HierL2'
% Revised: 12/10/13
%   -rewrote entirely, based on (new param) params.MODEL
%   -uses the matrix form of s (Ndims x Nmods)
% Repo'ed: 05/31/11
%   -this fxn is now used again
% Revised: 01/31/11
%   -This fxn is no longer in use---though it seems useful.  
%   -Now each column of a particular slice is in the same space *and* the
%   same modality.
% Cribbed: 01/10/11
%   from estError.m
%   by JGM
%-------------------------------------------------------------------------%


% set the neutral equal to the non-neutral (fix the rest below)
sN = s;

% convert non-neutral to neutral
switch params.MODEL
    case {'2Dinteg','1Dinteg'}
        switch params.NS
            case 'Hand-Position'
                propLogInd = strcmp(params.mods,'Joint-Angle');
                sN(:,propLogInd) = FK2link(s(:,propLogInd)',params,0);
            case 'Joint-Angle' 
                visLogInd = strcmp(params.mods,'Hand-Position');
                sN(:,visLogInd) = IK2link(s(:,visLogInd)',params,0);
            otherwise
                error('incorrect params.NS for this model');
        end
    case {'2Daddition','1Daddition'}
        switch params.NS
            case {'Hand-Position','Gaze-Angle'}
                propLogInd = strcmp(params.mods,'Joint-Angle');
                sN(:,propLogInd) = FK2link(s(:,propLogInd)',params,0);
            otherwise
                error('incorrect params.NS for addition models');
        end
    case 'HierL2'
        % do nothing---but NB that this expects the stimuli s fed into this
        % function to be purely in the neutral space
    otherwise
        error('unknown model for estGather')
end
        
    

end
