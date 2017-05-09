function yrdir = getdir(dirtype)
% getdir    Return a string with this computer's directory of type dirtype

%-------------------------------------------------------------------------%
% Revised: 07/28/16
%   -added 'tikz' dirs
% Revised: 06/28/16 (happy b'day TRM)
%   -renamed from getDataDir to getdir
%   -added argument; now you can ask for data dir OR code dir
% Revised: 04/18/16
%   -lowercased 'Themistocles'
% Created: 04/??/16
%   by JGM
%-------------------------------------------------------------------------%

[~,machine] = system('hostname');
switch dirtype
    case 'data'
        switch strtrim(machine)
            case {'CUPCAKE','Themistocles'}
                yrdir = 'C:\#DATA\';
            case {'PEPPERONI','ANCHOVY','MUSHROOM','ZAMFIR'}
                yrdir = 'C:\Users\makin\Data\';
            case 'domestica'
                yrdir = '~/data/';
            otherwise
                %%%error('unrecognized machine! -- jgm');
                fprintf('unrecognized machine! -- assuming domestica...\n');
                yrdir = '~/data/';
        end
    case 'code'
        switch strtrim(machine)
            case {'CUPCAKE','Themistocles'}
                yrdir = 'C:\#code\';
            case {'PEPPERONI','ANCHOVY','ZAMFIR','MUSHROOM'}
                yrdir = 'C:\Users\makin\code\';
			case 'domestica'
				yrdir = '~/code/';
            otherwise
                %%%error('unrecognized machine! -- jgm');
                fprintf('unrecognized machine! -- assuming domestica...\n');
                yrdir = '~/code/';
        end
    case 'tikz'
        switch strtrim(machine)
            case 'kobayashi-maru'
                yrdir = 'C:\Documents and Settings\makin\My Documents\#texs\tikzpics\';
            case {'CUPCAKE','Themistocles'}
                yrdir = 'C:\Users\makin\Documents\#texs\tikzpics\';
            case {'MUSHROOM','ANCHOVY','pepperoni','zamfir'}
                yrdir = 'C:\Users\makin\Documents\';
            case 'domestica'
                yrdir = '~/tikzpics/';
            otherwise
                %%% error('unrecognized machine! -- jgm');
                fprintf('unrecognized machine! -- assuming domestica...\n');
                yrdir = '~/tikzpics/';
        end
    otherwise
        error('unrecognized directory type -- jgm');
end