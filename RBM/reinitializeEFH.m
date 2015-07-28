function [vishid,hidbiases,visbiases,vishidinc,hidbiasinc,visbiasinc] =...
    reinitializeEFH(i_rbm,numsUnits,wts)
% reinitializeEFH   Initialize weights, biases, and increments of EFH

%-------------------------------------------------------------------------%
% Revised: ??/??/14
%   -changed allocation for GPU compatibility
% Created: ??/??/?? (very early)
%   by JGM
%-------------------------------------------------------------------------%


% init
Nvis = numsUnits(i_rbm);
Nhid = numsUnits(i_rbm+1);
[~,machine] = system('hostname');
machine = strtrim(machine);
if strcmp(machine,'domestica'),
    yrclass = 'gpuArray';
else
    yrclass = 'double';
end

% initialize symmetric weights and biases.
vishid      = 0.01*randn(Nvis, Nhid,yrclass);
hidbiases   = 0*ones(1,Nhid,yrclass);
% hidbiases   = -2*ones(1,Nhid,yrclass);
visbiases   = 0*ones(Nvis,1,yrclass);

if i_rbm > 1
    if Nhid==numsUnits(i_rbm-1)
        numRBMs     = length(numsUnits)-1;
        vishid      = wts{2*numRBMs-i_rbm+2}(1:end-1,:);% i.e. hidvis
        hidbiases   = wts{2*numRBMs-i_rbm+2}(end,:);    % i.e. the visbiases
        visbiases   = wts{i_rbm-1}(end,:)';             % i.e. the hidbiases
    end
end

vishidinc	= zeros(Nvis,Nhid,yrclass);
hidbiasinc  = zeros(1,Nhid,yrclass);
visbiasinc  = zeros(Nvis,1,yrclass);



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
