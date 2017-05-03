function [R,Q] = getDataMoCap(S,Q,testortrain)
% getDataMoCap
%  
% Get motion capture data for a EFH training.  S and Q can be empty, but S
% must have size (Ntraj x M x T), where M can (and should, to save memory)
% be 0.
%
% Standard use is to wrap this in an anonymous function with the params,
% and then pass to generateData.m

%-------------------------------------------------------------------------%
% Cribbed: 01/01/17
%   from generateData.m
%   by JGM
%-------------------------------------------------------------------------%

% "Ns"
Nexamples = size(S,1);
T = Q.T;
Ntraj = floor(Nexamples/T);

% hard-coded
TOPLOT = 0;
dataset = 'default';  % 'Jog';
test_frac = 0.2;

% load raw mocap data and naively drop frames (after I.S.)
switch dataset
    case 'default'
        load([getdir('data'),'MoCap/data'],'skel','Motion');
        for iSeq=1:length(Motion)
            Motion{iSeq} = cast(Motion{iSeq}(1:4:end,:),'like',S);
        end
    case 'Jog'
        load([getdir('data'),'MoCap/Jog1_M'],'skel','Walking');
        skel.type = 'mit';
        for iSeq=1:length(Walking)
            Motion{iSeq} = cast(Walking{iSeq}(1:4:end,:),'like',S);
        end
    otherwise
        error('whoops!');
end


% separate training from testing data
for iSeq = 1:length(Motion)
    Lseq = size(Motion{iSeq},1);
    Ntest = round(test_frac*Lseq);
    switch testortrain
        case 'train',   i0 = 1;             iF = Lseq-Ntest;
        case 'test',    i0 = Lseq-Ntest+1;  iF = Lseq;
        otherwise
            error('whoops!')
    end
    Motion{iSeq} = Motion{iSeq}(i0:iF,:);
end


% convert to EFH-like data (see help)
data = mocapBase2EFHdata(skel,Motion);

% plot?
if TOPLOT
    figure(300);
    subplot(3,1,1); plot(data.R(data.seqindex{1},:))
    subplot(3,1,2); plot(data.R(data.seqindex{2},:))
    subplot(3,1,3); plot(data.R(data.seqindex{3},:))
    
    % back to original format for plotting
    Motion = mocapEFHdata2Base(skel,data);
    figure(200); expPlayData(skel, Motion{1}, 1/30)
end

% break up data into Ntraj-long chunks
startindSet = [];
for iSeq = 1:length(data.seqindex)
    startindSet = [startindSet, data.seqindex{iSeq}(1:end-T+1)];
end
randinds = ceil(length(startindSet)*rand(Ntraj,1));
startinds = startindSet(randinds);
allinds = startinds + (0:T-1)';
R = data.R(allinds,:);


end