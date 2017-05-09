function getRemoteData(filename,serverloc,clientloc)
% USAGE:
%   getRemoteData(serverloc,clientloc)
%

%-------------------------------------------------------------------------%
% Created 03/13/17 (happy b'day, VMW)
%   by JGM
%-------------------------------------------------------------------------%

if ~exist([clientloc,filename],'file')
    if strcmp(computer,'GLNXA64')
        command = ['scp makin@7layerburrito.cin.ucsf.edu:',...
            serverloc,filename,' ',clientloc,'.'];
        unix(command);
    else
        command = ['pscp -scp -unsafe makin@7layerburrito:',...
            serverloc,filename,' ',clientloc,'.'];
        dos(command);
        %%% DOS programs always return status = 0, even when they work
    end
    if isempty(ls([clientloc,filename]))
        error('file failed to download');
    end
else
    fprintf('using local copy of data file\n');
end


end