% Version 1.000
%
% Code provided by Geoff Hinton and Ruslan Salakhutdinov
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with a
% note saying that the original programs are available from our web page.
% The programs and documents are distributed without any warranty, express
% or implied.  As the programs were written for research purposes only,
% they have not been tested to the degree that would be advisable in any
% important application.  All use of these programs is entirely at the
% user's own risk.

% This program trains Restricted Boltzmann Machine in which visible,
% binary, stochastic pixels are connected to hidden, binary, stochastic
% feature detectors using symmetrically weighted connections. Learning is
% done with 1-step Contrastive Divergence. The program assumes that the
% following variables are set externally:
%
%   maxepoch    -- maximum number of epochs 
%   numhid      -- number of hidden units 
%   batchdata   -- the data that is divided into batches 
%   (numcases numdims numbatches) 
%   restart     -- set to 1 if learning starts from beginning
%   LINEAR      -- linear vs. binary-sigmoidal hidden units
%-------------------------------------------------------------------------%
% Revised: 5/25/10
%   -made number of contrastive-divergence steps variable
%   -incorporated rbmhidlinear
% Revised: 5/24/10
%   -fixed formatting
% modified by JGM
%-------------------------------------------------------------------------%

% params
DISP = 0;
HIDFXN = params.typeUnits{i_rbm+1};
VISFXN = params.typeUnits{i_rbm};
alpha = 100;                                % for learning linear units
epsilonw = params.epsilonw;
epsilonvb = params.epsilonvb;
epsilonhb = params.epsilonhb; % *(alpha^strcmp(HIDFXN,'Bernoulli')); 
weightcost  = params.weightcost;
numtestbatches = params.numtestbatches;
finalmomentum = params.finalmomentum;
initialmomentum = params.initialmomentum;
tErravgOld = inf;
counter = 0;

% init
[numcases numdims numbatches] = size(batchdata);
if restart == 1,
    restart = 0;
    epoch = 1;

    % Initializing symmetric weights and biases.
    vishid      = 0.01*randn(numdims, numhid);
    hidbiases   = zeros(1,numhid);
    visbiases   = zeros(1,numdims);  % 3.21*ones(1,numdims);

    poshidmeans = zeros(numcases,numhid);
    neghidmeans = zeros(numcases,numhid);
    posprods    = zeros(numdims,numhid);
    negprods    = zeros(numdims,numhid);
    vishidinc   = zeros(numdims,numhid);
    hidbiasinc  = zeros(1,numhid);
    visbiasinc  = zeros(1,numdims);
    batchposhidprobs = zeros(numcases,numhid,numbatches);
end

% cycle through training data
indices = ceil(numhid*rand(rows*cols,1));
for epoch = epoch:maxepoch,
    fprintf(1,'epoch %d\r',epoch);
    errsum = 0;
    
    parfor batch = 1:numbatches,
        fprintf(1,'epoch %d batch %d\r',epoch,batch);       
        
        % positive phase
        posdata = batchdata(:,:,batch);
        poshidmeans = feedforward(posdata,vishid,hidbiases,HIDFXN,params);
        batchposhidprobs(:,:,batch) = poshidmeans;
        posprods = posdata'*poshidmeans;
        poshidact = sum(poshidmeans);
        posvisact = sum(posdata);        
    
        % negative phase
        [negvismeans, neghidmeans] = CDstepper(poshidmeans,vishid,...
            visbiases,hidbiases,HIDFXN,VISFXN,params);
        negprods  = negvismeans'*neghidmeans;
        neghidact = sum(neghidmeans);
        negvisact = sum(negvismeans);
        
        % for printing (only)
        err = sum(sum((posdata - negvismeans).^2))/numcases;
        errsum = err + errsum;

        % ramp (up) momentum
        momentum = (epoch>5)*finalmomentum + (epoch<=5)*initialmomentum;
        
        % weight/bias update
        vishidinc = momentum*vishidinc +...
            epsilonw*((posprods - negprods)/numcases - weightcost*vishid);
        visbiasinc = momentum*visbiasinc +...
            (epsilonvb/numcases)*(posvisact - negvisact);
        hidbiasinc = momentum*hidbiasinc +...
            (epsilonhb/numcases)*(poshidact - neghidact);
        vishid = vishid + vishidinc;
        visbiases = visbiases + visbiasinc;
        hidbiases = hidbiases + hidbiasinc;

    end
    erravg = errsum/numbatches;
    
    
    % test
    [D2 x2test] = DATAGEN(numtestbatches,params);
    % testdata = D2/max(max(max(D2)));
    testdata = D2;
    tErrSum = 0;
    for i_testbatch = 1:numtestbatches
        tdata = testdata(:,:,i_testbatch);
        htest = feedforward(tdata,vishid,hidbiases,HIDFXN,params);
        vtest = feedforward(htest,vishid',visbiases,VISFXN,params);
        tErr = sum(sum((tdata - vtest).^2 ))/numcases;
        tErrSum = tErrSum + tErr;
    end
    tErravg = tErrSum/numtestbatches;
    if tErravg > tErravgOld
        counter = counter + 1;
    else
        counter = 0;
    end
    tErravgOld = tErravg;
    
    if counter > params.counterMax
        break
    end
    
    % say error
    fprintf(1,'epoch %4i error %6.4e terror %6.4e \n',epoch,erravg,tErravg);
    save rbmwts vishid hidbiases visbiases params epoch

end
