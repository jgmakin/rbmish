function wts = train(batchdata,wts,params)
% the training loop for a DBN

%-------------------------------------------------------------------------%
% Cribbed 10/02/12
%   -from recalibrate.m
%   by JGM
%-------------------------------------------------------------------------%

%%%%%%%%%
% None of this will work with the current code (02/24/16)
%   -JGM
%%%%%%%%%


% init
iRBM = 1;
numRBMs = 1;
vishid = wts{iRBM}(1:end-1,:);
hidbiases = wts{iRBM}(end,:);
visbiases = wts{numRBMs*2-iRBM+1}(end,:)';
visDstrbs = params.typeUnits{iRBM};
visNums = params.numsUnits{iRBM};
hidDstrbs = params.typeUnits{iRBM+1};
hidNums = params.numsUnits{iRBM+1};
Nhid = sum(hidNums);
[Ncases,Nvis,Nbatches] = size(batchdata);
vishidinc = zeros(Nvis,Nhid);
hidbiasinc = zeros(1,Nhid);
visbiasinc = zeros(1,Nvis);

% learning parameters
lrnrt = 10e-5; % 2e-5;
epsilonw =  lrnrt;
epsilonvb = lrnrt;
epsilonhb = lrnrt; 
momentum = params.initialmomentum;
weightcost = params.weightcost;

% train for one epoch
for batch = 1:Nbatches,
    fprintf('.');
    
    % positive phase
    posdata = batchdata(:,:,batch);
    poshidmeans = invParamMap(posdata,vishid,hidbiases,hidDstrbs,hidNums,params);
    poshidstates = sampleT(poshidmeans,hidDstrbs,hidNums,params);
    posprods = posdata'*poshidstates;
    poshidact = sum(poshidstates);
    posvisact = sum(posdata);
    
    % negative phase
    [negvisstates, neghidstates] = CDstepper(poshidstates,vishid,...
        visbiases,hidbiases,hidDstrbs,visDstrbs,hidNums,visNums,...
        params.Ncdsteps,params);
    negprods  = negvisstates'*neghidstates;
    neghidact = sum(neghidstates);
    negvisact = sum(negvisstates);
    
    % weight/bias update
    vishidinc = momentum*vishidinc +...
        epsilonw*((posprods - negprods)/Ncases - weightcost*vishid);
    visbiasinc = momentum*visbiasinc +...
        (epsilonvb/Ncases)*(posvisact - negvisact);
    hidbiasinc = momentum*hidbiasinc +...
        (epsilonhb/Ncases)*(poshidact - neghidact);
    vishid = vishid + vishidinc;
    visbiases = visbiases + visbiasinc;
    hidbiases = hidbiases + hidbiasinc;
    
end
% fprintf('\n');

wts{iRBM} = [vishid; hidbiases];
wts{numRBMs*2-iRBM+1} = [vishid'; visbiases'];


end