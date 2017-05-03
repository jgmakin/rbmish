% adaptinit.m
% set the path etc. for the adaptation code

%-------------------------------------------------------------------------%
% Revised: 09/27/12
%   -added check for hostname
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

[blank, name] = system('hostname');

switch strtrim(name)
    case 'kobayashi-maru'
        codedir = 'C:\#code';
        set(0,'DefaultFigurePosition',[1886 271 560 420]); clc;
    case 'CUPCAKE'
        set(0,'DefaultFigurePosition',[563 784 560 420]); clc;
    case {'mushroom','keck-phaser1','pepperoni','zamfir'}
        codedir = 'C:/Users/makin/code';
        set(0,'DefaultFigurePosition',[445 272 560 420]); clc;
    otherwise
        error('unknown host -- jgm');
end

path(path,[codedir,'/rbm']);
path(path,[codedir,'/rbm/scratch']);
path(path,[codedir,'/rbm/parallel']);
path(path,[codedir,'/rbm/retired']);
path(path,[codedir,'/rbm/results']);
path(path,[codedir,'/rbm/tuningcurves']);
path(path,[codedir,'/robotics']);
path([codedir,'/utils'],path);      % this directory contains your tex.m
path(path,[codedir,'/tools']);






