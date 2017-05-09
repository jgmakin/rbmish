function [R,Q] = getDataPolyphonicMusic(S,Q,database,testortrain)
% getDataPolyphonicMusic
%
% USAGE:
%   R = getDataPolyphonicMusic(S,Q,params.database,params.testortrain);y

%-------------------------------------------------------------------------%
% Created: 12/??/16
%   by JGM
%-------------------------------------------------------------------------%

% "Ns"
[Nseq,~,T] = size(S);

% which polyphonic-music data set?
temp = load([getdir('data'),'musicdata/',database,'.mat'],testortrain);
data = temp.(testortrain);

% break up data into Nseq-long chunks
Nsongs = length(data);
allLengths = arrayfun(@(ii)(size(data{ii},2)),1:Nsongs);
seqindices = mat2cell(1:sum(allLengths),1,allLengths);
startindSet = [];
for iSeq = 1:Nsongs
    startindSet = [startindSet, seqindices{iSeq}(1:max(end-T+1,0))];
end
randinds = ceil(length(startindSet)*rand(Nseq,1));
startinds = startindSet(randinds);
allinds = startinds + (0:T-1)';
alldata = cat(2,data{:});
note0 = 22;     % none of these pieces uses a note below this
R = alldata(note0:end,allinds)';


end