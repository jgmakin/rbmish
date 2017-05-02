% Version 1.000
%
% Code provided by Ruslan Salakhutdinov and Geoff Hinton
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

% This program fine-tunes an autoencoder with backpropagation. Weights of
% the autoencoder are going to be saved in mnist_weights.mat and trainig
% and test reconstruction errors in mnist_error.mat You can also set
% maxepoch, default value is 200 as in our paper.
%-------------------------------------------------------------------------%
% Revised: 08/03/10
%   -added params as an argument to uppass
%   -added params as an arg to CG_func (for the same reason)
% Revised: 5/25/10
%   -fixed formatting
%   -modified to accomodate variable number of layers
% modified by JGM
%-------------------------------------------------------------------------%

%init
maxepoch = 200;
max_iter = 3;                           % number of linesearches
numMinibatches = 10;                    % to be combined in a big batch
numClasses = 10;                        % there are 10 Arabic numerals
typeUnits = params.typeUnits;
numsUnits = params.numsUnits;

makebatches;
[Ncases,Nvis,Nbatches] = size(batchdata);
N = Ncases;
fprintf(1,'\nTraining discriminative model on MNIST by minimizing ');
fprintf(1,'cross-entropy error. \n');
fprintf(1,'\n%i batches of %i cases each. \n',Nbatches,Ncases);

% load weights (etc.)
clear wts
load ClassifWtsFile
wts{end+1} = 0.1*randn(size(wts{end},2)+1, numClasses);

% nodes/layer of the classifier network
unitvec = arrayfun(@(i)(sum(numsUnits{i})),1:length(numsUnits));
Dim = [unitvec numClasses]';

% init errors
test_err = zeros(maxepoch,1);
test_crerr = zeros(maxepoch,1);
train_err = zeros(maxepoch,1);
train_crerr = zeros(maxepoch,1);

% backprop
for epoch = 1:maxepoch

    % COMPUTE TRAINING MISCLASSIFICATION ERROR
    [train_err(epoch) train_crerr(epoch)] =...
        uppass(batchdata,batchtargets,wts,params);
    
    % COMPUTE TEST MISCLASSIFICATION ERROR
    [test_err(epoch) test_crerr(epoch)] =...
        uppass(testbatchdata,testbatchtargets,wts,params);
        
    % print
    [Ncases,Nvis,Nbatches] = size(batchdata);    
    [Ntestcases,Ntestvis,Ntestbatches] = size(testbatchdata);
    fprintf(1,'Before epoch %d Train # misclassified: ',epoch);
    fprintf(1,'%d (from %d). Test # misclassified: ',train_err(epoch),...
        Ncases*Nbatches);
    fprintf(1,'%d (from %d) \t \t \n',test_err(epoch),...
        Ntestcases*Ntestbatches);

    % conjugate gradient method
    i_bigbatch = 0;
    data = zeros(Ncases*numMinibatches,Nvis);        % malloc
    targets = zeros(Ncases*numMinibatches,numClasses);   % malloc
    for batch = 1:Nbatches/numMinibatches
        fprintf(1,'epoch %d batch %d\r',epoch,batch);
        i_bigbatch = i_bigbatch + 1;
            
        % COMBINE 10 MINIBATCHES INTO 1 LARGER MINIBATCH
        offset = 0;
        for kk = 1:numMinibatches
            data(offset+1:offset+Ncases,:) =...
                batchdata(:,:,(i_bigbatch-1)*numMinibatches + kk);
            targets(offset+1:offset+Ncases,:) =...
                batchtargets(:,:,(i_bigbatch-1)*numMinibatches + kk);
            offset = offset + Ncases;
        end
        
        % PERFORM CONJUGATE GRADIENT WITH max_iter LINESEARCHES
        % First update top-level weights holding other weights fixed.
        if epoch < 6
            
            % push data through the network (but not through output)
            probs = data;
            for layer = 1:length(wts)-1
                probs = invParamMap(probs,wts{layer}(1:end-1,:),...
                    wts{layer}(end,:),typeUnits{layer+1},...
                    numsUnits{layer+1},params);
            end
            
            % update *class weights* only, with conjugate-gradient
            [wtsVec, fX] = minimize(wts{end}(:),'CG_func',max_iter,...
                Dim(end-1:end),params,probs,targets,'CLASSIFY');

            % unvectorize class weights
            wts{end} = reshape(w_classVec,Dim(end-1)+1,Dim(end));
        else
            
            % vectorize all weights
            wtsVec = [];
            for i = 1:length(wts)
                wtsVec = [wtsVec; wts{i}(:)];
            end
                       
            % update all weights with conjugate-gradient method
            [wtsVec, fX] = minimize(wtsVec,'CG_func',max_iter,...
                Dim,params,data,targets,'CLASSIFY');
            
            % unvectorize all (updated) weights
            offset = 0;
            for i = 1:length(wts)
                wts{i} = reshape(wtsVec(offset+1:offset+(Dim(i)+1)*...
                    Dim(i+1)), Dim(i)+1, Dim(i+1));
                offset = offset + (Dim(i)+1)*Dim(i+1);
            end            
        end
        % END OF CONJUGATE GRADIENT WITH 3 LINESEARCHES 

    end

    save mnistclassify_weights wts % w_class is wts{end}
    save mnistclassify_error test_err test_crerr train_err train_crerr;

end



