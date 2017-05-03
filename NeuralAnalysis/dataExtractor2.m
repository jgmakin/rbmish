function [X,y,angle] = dataExtractor2(day,trgnum,monkey)
% DATAEXTRACTOR2    Data extractor for HHS's original project
%
%

%-------------------------------------------------------------------------%
% Modifed: 03/22/11
%   -made day and angle# variable, to be set outside this script
% Created: 03/15/11
%   by JGM
%-------------------------------------------------------------------------%

load HHS_DataForJGM

% separate monkeys
i_e = 1;
i_d = 1;
for i = 1:length(DATA_STRUCT)
    if strncmp(DATA_STRUCT(i).date,'E',1)
        ezra(i_e) = DATA_STRUCT(i);
        i_e = i_e + 1;
        
    elseif strncmp(DATA_STRUCT(i).date,'D',1)
        dmitri(i_d) = DATA_STRUCT(i);
        i_d = i_d + 1;
    else 
        error('unrecognized data label! -- jgm\n');
        
    end
end

switch monkey
    case 'E'
        maxdays = length(ezra);
        monkeydata = ezra;
    case 'D'
        maxdays = length(dmitri);
        monkeydata = dmitri;
    otherwise
        error('unknown monkey! -- jgm\n');
end

% which day, which target?
if day <= maxdays
    data = monkeydata(day);    
    
    % X = behavioral data, U = neural data
    y = data.initangs(data.exp_idx_bytrg{trgnum});
    X = data.ratemat(data.exp_idx_bytrg{trgnum},:);
    
    % check
    angles = data.trgangs(data.exp_idx_bytrg{trgnum});
    if sum(angles == angles(1)) < length(angles)
        error('indices are confuzored -- jgm');
    else
        angle = angles(1);
    end
    
else
    fprintf('reached maximum number of days\n\n');
    angle = inf; X = inf; y = inf;
end


end