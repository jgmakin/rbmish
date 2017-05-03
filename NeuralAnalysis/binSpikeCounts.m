function [R,Xout,endinds] = binSpikeCounts(S,UnitSpikes,params)
% Bin spike counts, averaging the underlying dynamics

%-------------------------------------------------------------------------%
% Revised: 05/20/13
%   -more efficient data storage (outputs inds rather than XPast, XFuture)
% Cribbed: 05/16/13
%   -from KF4HHS
%   by JGM
%-------------------------------------------------------------------------%

% init params
Nneurons = length(UnitSpikes);
Ntrials = length(S);
Nstates = params.Nstates;
Ndims = params.Ndims;
dt = params.dt;
m = params.m;                               % 16 => 60 Hz (66.7 ms bins)
binsize = m*dt;                             %

% initialize
endinds = zeros(Ntrials+1,1,'like',S(1).t);
ind = 0;
for iTrial = 1:Ntrials
    if strcmp(params.BINMETHOD,'slidingwindow')
        Nbins = length(S(iTrial).t);
    else
        edges = S(iTrial).t(1):binsize:S(iTrial).t(end);
        Nbins = length(edges) - 1;
    end
    ind = ind+Nbins;
    endinds(iTrial+1) = ind;
end


% malloc
R = zeros(ind,Nneurons,'like',S(1).t);
Xout = zeros(ind,Nstates,'like',S(1).t);

% loop through trials, collecting up binned spike rates
for iTrial = 1:Ntrials
    
    
    if strcmp(params.BINMETHOD,'slidingwindow')
        
        % set the edges of the bins
        EdgeMat = getEdgeMatrix(S(iTrial).t(1),S(iTrial).t(end),dt,m);
        Nbins = size(EdgeMat,2)-1;
        
        % malloc
        Spks = zeros(m,Nbins,Nneurons);
        for iSlide = 1:m
            Spks(iSlide,:,:) = histcAllNeurons(UnitSpikes,EdgeMat(iSlide,:));
        end
        
        % convert to spike *rate* and store
        binsizes = repmat(diff(EdgeMat,[],2),[1,1,Nneurons]);
        theseRates = reshape(Spks./binsizes,m*Nbins,Nneurons);
    else
        
        % bin across all neurons
        Spks = histcAllNeurons(UnitSpikes,...
            S(iTrial).t(1):binsize:S(iTrial).t(end));
        
        % convert to spike *rate*
        %%%theseRates = Spks/binsize;
        theseRates = Spks;
        %%% why bother converting to rates?
        Nbins = length(S(iTrial).t(1):binsize:S(iTrial).t(end))-1;
    end
    
    % store the (averaged/downsampled) dynamics during this trial
    if isfield(S,'pos')
        switch Nstates/Ndims %%% not very elegant
            case 1
                thisX = [S(iTrial).pos];
            case 2
                thisX = [S(iTrial).pos S(iTrial).vel];
            case 3
                thisX = [S(iTrial).pos S(iTrial).vel S(iTrial).acc];
        end
    else
        thisX = S(iTrial).X;
    end
        
        
    switch params.BINMETHOD
        case 'downsampled'
            thisXout = thisX(m:m:(Nbins*m),:);
        case 'slidingwindow'
            thisXout = thisX;
        otherwise
            thisXout = squeeze(mean(reshape(thisX(1:(Nbins*m),:)',...
                [Nstates,m,Nbins]),2))';
    end
    Xout((endinds(iTrial)+1):(endinds(iTrial+1)),:) = thisXout;
         
    % store the outputs
    R((endinds(iTrial)+1):(endinds(iTrial+1)),:) =...
        cast(theseRates(1:size(thisXout,1),:),'like',Ndims);
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function Spks = histcAllNeurons(UnitSpikes,edges)
% just applies histc to all the neurons in the struct UnitSpikes

Spks = arrayfun(@(ii)(histc(UnitSpikes(ii).t(:)',edges(:)')),...
    1:length(UnitSpikes),'UniformOutput',false);
Spks = cat(1,Spks{:})';
Spks = Spks(1:end-1,:); % trim off last bin

end
%-------------------------------------------------------------------------%


