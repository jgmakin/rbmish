function [X,Q] = getLatentsCoherences(Nexamples,yrclass,params,varargin)
% getLatentsCoherences
%
% USAGES:
%   [X,Q] = getLatentsCoherences(Nexamples,yrclass,params)
%
%   [X,Q] = getLatentsCoherences(Nexamples,yrclass,params,T)

%-------------------------------------------------------------------------%
% Cribbed: 01/02/17
%   from generateData.m
%   by JGM
%-------------------------------------------------------------------------%

% load and clip coherences
trainsuffix = params.trainsuffix(ceil(size(params.trainsuffix,1)*rand),:);
Q.coh = mean(loadAndClipCoherences(params.subj,trainsuffix,params),2);
%%%%%%%
% Nsamples = length(Q.coh);
% N = 1*60*params.Nsperwindow;
% while length(Q.coh) <= 1*Nbatches
% 	trainsuffix = params.trainsuffix(...
%   ceil(size(params.trainsuffix,1)*rand),:);
%   Q.coh = mean(loadAndClipCoherences(params.subj,trainsuffix,params),2);
% end
% iii = (1:N:floor((Nsamples-N)/N)*N)' + (1:N)-1;
% Q.coh = exp(mean(log(Q.coh(iii)),2));
%%%%%%%


% Ns
Nvis = params.numsUnits{1}(end); %%% something of a hack
if strcmp(params.typeUnits{1}(end),'GammaFixedScale')
    Ncoh  = Nvis;
else
    Ncoh = Nvis/2;              % the other half are log(coh)
end
L = length(Q.coh);
%%%winsize = 12;
%%%winstep = 6;
% b/c the coherence were themselves computed in 10 second
% windows slid by 5 seconds, 12 of them gives 1 min.
% [[CHECK THESE NUMBERS WITH LK!!]]
Nminperinterval = 0.25; %%%10;           % interval for coh. var.: 10min
Nwinperinterval = round(Nminperinterval*60/params.Nsperwinstep);
Nwinperintervalstep = round(Nwinperinterval/2);


% compute the variance in a sliding window
Nchunks = floor((L - Nwinperinterval)/Nwinperintervalstep);
MaxInd = Nchunks*Nwinperintervalstep;
inds = (1:Nwinperintervalstep:MaxInd) + ((1:Nwinperinterval)-1)';
CohVar = var(Q.coh(inds));

% make explicit assumption that coherence is const. w/i window
CohVar = repmat(CohVar,[Nwinperintervalstep,1]);
try
    CohVar = [CohVar(:);repmat(CohVar(end),[L-MaxInd,1])];
catch me
    fprintf('making up coherence variances...\n');
    CohVar = rand(size(Q.coh));
end
%plotyy(1:100,coherenceTraj(1:100),1:100,CohVar(1:100));



%%%%
% consider replacing CohVar with something based on cdf....
%%%%
% stimuli and coherence indices
iSequencelength = find(strcmp(varargin,'sequencelength'));
if ~isempty(iSequencelength)
    
    % "Ns"
    T = varargin{iSequencelength + 1};
    Ntraj = floor(Nexamples/T);
    
    % construct inds such that the data are *consecutive*
    Nconsecutive = Ncoh*T;
    inds = ceil((L-Nconsecutive+1)*rand(Ntraj,1));
    inds = inds + (0:(Nconsecutive-1));
    inds = reshape(inds,[Ntraj,Ncoh,T]);
    
    % the "stim" is whether we're in a low or high var state
    X = double(mean(CohVar(inds),2) > mean(CohVar));
    X = longdata(X);
    inds = longdata(inds);
    
    % 
    Q.restarts  = 1:T:Nexamples;
    Q.T         = T;
    
else
    % consecutive "trials" are uncorrelated
    inds = ceil((L-Ncoh+1)*rand(Nexamples,1));
    inds = inds + (0:(Ncoh-1));
    
    % the "stim" is whether we're in a low or high var state
    X = double(mean(CohVar(inds),2) > mean(CohVar));
end

if strcmp(yrclass,'gpuArray'), X = gpuArray(X); end

% store for use in generating R
Q.inds = inds;




end