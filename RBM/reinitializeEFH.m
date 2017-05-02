function [vishid,hidbiases,visbiases,vishidinc,hidbiasinc,visbiasinc] =...
    reinitializeEFH(iRBM,numsUnits,typeUnits,wts,datatype,dataclass)
% reinitializeEFH   Initialize weights, biases, and increments of EFH

%-------------------------------------------------------------------------%
% Revised: 12/28/16
%   -eliminated inputs and outputs concerning recurrent indices, types,
%   nums.
% Revised: 03/23/16
%   -incorporated EFHinitNs.m (which is now retired)
% Revised: ??/??/14
%   -changed allocation for GPU compatibility
% Created: ??/??/?? (very early)
%   by JGM
%-------------------------------------------------------------------------%


% init
visDstrbs = typeUnits{iRBM};
visNums = numsUnits{iRBM};
Nvis = sum(visNums);
Nhid = sum(numsUnits{iRBM+1});

% (re)initialize the wts/biases and wt/bias velociities
if iRBM==1
    switch datatype
        case {'ErlangMixtureToy','ECcoherences'}
            fprintf('\n\n\nSpecial weight initialization for ');
            fprintf('ErlangMixture datatype!!\n\n\n');
            %vishid(1:Nvis/2,:) = repmat(0.005*randn(1,Nhid,yrclass),[Nvis/2,1]);
            %vishid((Nvis/2+1):Nvis,:) = repmat(0.005*randn(1,Nhid,yrclass),[Nvis/2,1]);
            % you could tie the weights, but it's not apparently necessary
            vishid      = 0.005*randn(Nvis,Nhid,dataclass);
            hidbiases   = 0*ones(1,Nhid,dataclass);
            visbiases   = 0*ones(Nvis,1,dataclass);
            if strcmp(visDstrbs{1}(end),'Erlang')||strcmp(visDstrbs{1}(end),'Gamma')
                visbiases((end/2+1):end) = -4;
            end
        otherwise
            % initialize symmetric weights and biases.
            vishid      = 0.005*randn(Nvis, Nhid,dataclass); %%% 0.005 for I.S., I think
            hidbiases   = 0*ones(1,Nhid,dataclass);
            visbiases   = 0*ones(Nvis,1,dataclass);
    end
else
    if Nhid==sum(numsUnits{iRBM-1})
        numRBMs     = length(numsUnits)-1;
        vishid      = wts{2*numRBMs-iRBM+2}(1:end-1,:);% i.e. hidvis
        hidbiases   = wts{2*numRBMs-iRBM+2}(end,:);    % i.e. the visbiases
        visbiases   = wts{iRBM-1}(end,:)';             % i.e. the hidbiases
    else
        % initialize symmetric weights and biases.
        vishid      = 0.005*randn(Nvis, Nhid,dataclass); %%% 0.005 for I.S., I think
        hidbiases   = 0*ones(1,Nhid,dataclass);
        visbiases   = 0*ones(Nvis,1,dataclass);
    end
end
vishidinc	= zeros(Nvis,Nhid,dataclass);
hidbiasinc  = zeros(1,Nhid,dataclass);
visbiasinc  = zeros(Nvis,1,dataclass);



%%%
% load('dynamical\finalwts\wts1DrEFHManyXprmts.mat','Allwts');
% vishid = Allwts{1}{1}(1:end-1,:);
% hidbiases = Allwts{1}{1}(end,:);
% visbiases = Allwts{1}{2}(end,:)';
%%%

%%% when did you ever use these??
% poshidmeans     = zeros(numcases,numhid);
% neghidstates    = zeros(numcases,numhid);
% posprods        = zeros(numdims,numhid);
% negprods        = zeros(numdims,numhid);


end
