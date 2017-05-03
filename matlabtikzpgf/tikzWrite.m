function tikzWrite(outtxt,tikzfilename,varargin)
% tikzWrite     Write out a tikz file
%
%   tikzWrite(outtxt,tikzfilename)
%
% Write the string outtxt into the file tikzfilename.

%-------------------------------------------------------------------------%
% Cribbed: 02/28/16
%   -by JGM
%-------------------------------------------------------------------------%

% wrap in tikzpicture
outtxt = ['\begin{tikzpicture}',varargin{:},sprintf('\n'),outtxt];
outtxt = [outtxt,sprintf('\n'),'\end{tikzpicture}%'];

% write to file
[blank, name] = system('hostname');
switch strtrim(name)
    case 'kobayashi-maru'
        yrtikzdir = 'C:\Documents and Settings\makin\My Documents\#texs\tikzpics\';
    case {'CUPCAKE','Themistocles'}
        yrtikzdir = 'C:\Users\makin\Documents\#texs\tikzpics\';
    case {'MUSHROOM','keck-phaser1','pepperoni','zamfir','ANCHOVY'}
        yrtikzdir = 'C:\Users\makin\Documents\';
    case 'domestica'
        yrtikzdir = '~/tikzpics/';
    otherwise
        error('unknown host -- jgm');
end
outfile = [yrtikzdir,tikzfilename,'.tex'];
fid = fopen(outfile,'wt+');
fwrite(fid,outtxt);
fclose(fid);

end
