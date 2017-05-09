function [R,Q] = getDataECoG(S,Q,machine)

% getDataECoG   Raw ECoG data for EFH training/testing
% 
% USAGE:
%   R = getDataECoG(S,Q,params)

%-------------------------------------------------------------------------%
% Cribbed: 12/31/16
%   by JGM
%   from generateData.m
%-------------------------------------------------------------------------%



% init
subj = 'EC108';
filesuffix = '_2fd7f1ac';
iFile = ceil(100*rand); %%% why not 144? or nearby?

filepredir = [subj,filesuffix,'/mat_allChannels/'];
filename = [subj,filesuffix,'_all_',sprintf('%.3d',iFile),'.mat'];
serverloc = ['/sabes_data3/kderosier/Chang_lab_data/',subj,'/'];
clientloc = [' ',getdir('data'),'SUBNETS/',subj,filesuffix,...
    '_KHPDpreprocessed/mat_allChannels/'];
fprintf('...attempting to load %s...\n\n',filename);
switch machine
    case 'domestica'
        if ~exist([clientloc(2:end),filename],'file')
            command = ['scp makin@7layerburrito.cin.ucsf.edu:',...
                serverloc,filepredir,filename,clientloc,'.'];
            status = unix(command);
        end
        load([clientloc(2:end),filename],'data')
        if ~exist('data','var')
            load([clientloc(2:end),filename],'smallData')
            data = smallData;
            clear smallData;
        end
        
    otherwise
        if exist(['Z:/',filepredir,filename],'file')
            load(['Z:/',filepredir,filename]);
        else
            command = ['pscp -scp makin@7layerburrito:',...
                serverloc,filepredir,filename,clientloc,'.'];
            status = dos(command);
            load([clientloc(2:end),filename],'data')
            if ~exist('data','var')
                load([clientloc(2:end),filename],'smallData')
                data = smallData;
                clear smallData;
            end
        end
end
%%%%%
% and now delete the file
%%%%%

indF = size(data,2);
if size(S,1) > indF
    error('whoops -- too many data requested! -- jgm');
end


% Don't bother to get Ntraj different trajectories, and then apply
% longdata--we've got plenty of data.  So just make one long trajectory of
% Nexamples consecutive samples (yes, this stores the data appropriately).
[Ntraj,~,T] = size(S);
Nexamples = Ntraj*T;
data = data(:,(0:Nexamples-1) + ceil((indF-Nexamples)*rand))';


% filter and normalize variance
fs = 512;   %%% hard-coded
fLow = 10;  %%% 10
fHi = 50;   %%% 80
data = eegfilt(data',fs,fLow,fHi,0,3*fix(fs/fLow),0,'fir1')';
[aminds,hpinds] = getAMHPinds(subj);
data = data(:,[aminds,hpinds]);
try
    % de-mean and whiten
    % Sigma = cov(data);
    % Sigma = Sigma*Sigma'; %%% for numerical stability
    % data = (data - mean(data))/chol(Sigma);
    R = (data - mean(data))./std(data);
catch ME
    keyboard
end
%%% data = data*3;
%%% *Some* of the variance is from the dynamics....


end