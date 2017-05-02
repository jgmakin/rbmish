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
% the autoencoder are going to be saved in mnist_weights.mat and training
% and test reconstruction errors in mnist_error.mat You can also set
% maxepoch; default value is 200 as in our paper.
%-------------------------------------------------------------------------%
% Revised: 08/03/10
%   -added params as an argument to updown
%   -added params as an arg to CG_func (for the same reason)
% Revised: 5/25/10
%   -modified to accomodate variable number of layers
% Revised: 5/24/10
%   -fixed formatting
% modified by JGM
%-------------------------------------------------------------------------%

% init
maxepoch = params.BPmaxepoch;
max_iter = params.max_iter;                 % number of linesearches
numMinibatches = params.numMinibatches;     % to be combined in a big batch
DISP = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%
% makebatches;
%%%%%%%%%%%%%%%%%%%%%%%%%
batchdata = batchdata0;
testbatchdata = testbatchdata0;


[Ncases,Nvis,Nbatches] = size(batchdata);
N = Ncases;
fprintf(1,'\nFine-tuning deep autoencoder by minimizing cross-entropy');
fprintf(1,' error. \n%i batches of %i cases each. \n',Nbatches,Ncases);
    
% load weights (etc.)
while 1
    r = input('clear weights? y/n\n','s');
    if strcmp(r,'y')
        clear wts
        load EncoderWtsFile
    elseif strcmp(r,'n')
        break;
    else
    end
end

numRBMs = length(numsUnits)-1;

% nodes/layer of the "unfolded" network
unitvec = arrayfun(@(i)(sum(numsUnits{i})),1:length(numsUnits));
Dim = unitvec([1:end,end-1:-1:1])';

% init errors
test_err = zeros(maxepoch,1);
train_err = zeros(maxepoch,1);

% backprop
for epoch = 1:maxepoch

    % COMPUTE TRAINING RECONSTRUCTION ERROR
    [train_err(epoch),dataout] = updownfast(batchdata, wts, params); 

    
    if DISP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        numDispCases = 15;
        griddisp(batchdata(:,:,end),dataout,numDispCases,params);
        
        
        % display
%         fprintf(1,'Displaying in figure 1: Top row - real data, ');
%         fprintf(1,'Bottom row -- reconstructions \n');
%         data = batchdata(1:numDispCases,:,end);     % final batch
%         Mdisp = [data dataout]';
%         Mdisp = reshape(Mdisp,size(data,2),numDispCases*2);
%         if epoch==1
%             close all
%             figure('Position',[100,600,1000,200]);
%         else
%             figure(1)
%         end
%         mnistdisp(Mdisp);
%         drawnow;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % COMPUTE TEST RECONSTRUCTION ERROR
    [test_err(epoch),dataout] = updownfast(testbatchdata, wts, params);  
    
    % print
    fprintf(1,'Before epoch %d Train squared error: ',epoch);
    fprintf(1,'%6.3f Test squared error: ',train_err(epoch));
    fprintf(1,'%6.3f \t \t \n',test_err(epoch));

    
    % conjugate gradient
    i_bigbatch = 0;
    data = zeros(Ncases*numMinibatches,Nvis);  % malloc
    for batch = 1:Nbatches/numMinibatches
        fprintf(1,'epoch %d batch %d\r',epoch,batch);
        i_bigbatch = i_bigbatch + 1;
        
        % COMBINE numMinibatches MINIBATCHES INTO 1 LARGER MINIBATCH
        offset = 0;
        for kk = 1:numMinibatches
            data(offset+1:offset+Ncases,:) =...
                batchdata(:,:,(i_bigbatch-1)*numMinibatches + kk);
            offset = offset + Ncases;
        end         
%         data = reshape(...
%             batchdata(:,:,i_bigbatch:i_bigbatch+numMinibatches-1),...
%             numcases*numMinibatches,numdims);
%         i_bigbatch = i_bigbatch + numMinibatches;
         

        % "minimize" (apparently) requires its "starting point"---the 
        %   weights---to be a vector 
        wtsVec = [];
        for i = 1:length(wts)
            wtsVec = [wtsVec; wts{i}(:)];
        end
        
        % run conjugate-gradient method to get better weights
        [wtsVec, fX] = minimize(wtsVec,'CG_func',max_iter,...
            Dim,params,data,'ENCODE');
        
        % put the (updated) weights back into matrices  
        offset = 0;
        for i = 1:length(wts)
            wts{i} = reshape(wtsVec(offset+1:offset+(Dim(i)+1)*...
                Dim(i+1)), Dim(i)+1, Dim(i+1));
            offset = offset + (Dim(i)+1)*Dim(i+1);
        end
    end

    save caca_weights wts
    save caca_error test_err train_err;

end



