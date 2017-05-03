% recalibration3: test recalibration for different initial shifts

%-------------------------------------------------------------------------%
% Created: 09/??/12
%   by JGM
%-------------------------------------------------------------------------%

clear all;

% load wts and params, fix params
% load ../results/numhidswts/Std050
load ../results/StdAllGains.mat
% load ../results/numhidswts/Std050
params.Ncases = 5; % 100;
params.smpls = 15;

% init
% gainvecs = [12 12; 12 18; 15 15; 18 12; 18 18];
gainvecs = [9 9]; % [03 03; 03 15; 09 09; 15 03; 15 15];
% gainvecs = [03 15; 09 09; 15 15];
ext = 'killme'; % 'quick'; % 'slow';
nSubjects = 8; % 100;
nBatches = 120;
% stddevs = 5;


% set how much prop is shifted from vis (in prop space)
shft = [0.05; 0.05]; % (radians)
for kk = 1:size(gainvecs,1)
    
    % run the recalibration inner loop
    params.gmin = gainvecs(kk,:);
    params.gmax = gainvecs(jj,:);
    [IntegL0,DoBF,stVWDR,stWDR,stEMP] =....
        recalibrationCorePP(nSubjects,nBatches,shft,wts,params);
    
    % save
    filename = ['adpt',num2str(params.gmin,'%02.f'),ext];
    save(filename,'IntegL0','DoBF','stVWDR','stWDR','stEMP','params');
    
end