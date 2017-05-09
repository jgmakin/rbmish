function clrName = getColorName(name)
% getColorName  Get color name from associated "tag"
% 
% USAGE:
%   clrName = getColorName(tags)
%
% The clrName should match some color name defined in rbmish.sty.  It is
% then used in some tikz plot.

%-------------------------------------------------------------------------%
% Revised: 03/30/15 (happy b'day, DMM)
%   -added case 'Hand-Velocity'
% Created: 10/22/14
%   by JGM
%-------------------------------------------------------------------------%

setColors

switch name
    case 'Hand-Position'
        clrName = 'vclr';
    case 'Joint-Angle'
        clrName = 'pclr';
    case 'Gaze-Angle'
        clrName = 'eyeclr';
    case 'Efference-Copy'
        clrName = 'efcpclr';
    case {'rEFH','EFH'}
        clrName = 'rbmclr';
    case 'opt'
        clrName = 'optclr';
    case 'EM'
        clrName = 'EMclr';
    case 'EM (best)'
        clrName = 'EMBESTclr';
    case 'obs'
        clrName = 'OBSclr';
    case {'obs (no u)','obs 1st-ord'}
        clrName = 'OBSclr';
    case 'prior'
        clrName = 'OBSclr';
    case {'Hand-Velocity','Angular-Velocity'}
        clrName = 'velclr';
    otherwise
        fprintf('warning: unknown estimator name -- jgm\n');        
        clrName = 'OBSclr';
end
