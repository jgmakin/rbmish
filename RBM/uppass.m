function [train_err,train_crerr] =uppass(batchdata,batchtargets,wts,params)
% The counterpart for classification of the encoder's updown function for
% fine-tuning: propagate probabilities up through the network, keeping
% track of error...

%-------------------------------------------------------------------------%
% Revised: 08/03/10
%   -fixed the FXNs to refer to the right nodes
% Revised: 06/03/10
%   -changed feedforward FXN to strings
% Adapted: 05/25/10
%   from the Salakhutdinov/Hinton code
%   by JGM
%-------------------------------------------------------------------------%

% init
err_cr = 0;
counter = 0;
[Ncases,~,Nbatches] = size(batchdata);

% loop through batches
for batch = 1:Nbatches
    data = batchdata(:,:,batch);
    target = batchtargets(:,:,batch);

    % push up through layers (up pass only)
    stdParams = data;
    for layer = 1:length(wts)
        hidDstrbs = params.typeUnits{layer+1};
        hidNums = params.numsUnits{layer+1};
        stdParams = invParamMap(stdParams,wts{layer}(1:end-1,:),...
            wts{layer}(end,:),hidDstrbs,hidNums,params);
    end
    dataout = stdParams;

    % count how often the max. output is at the same digit as the target
    [~, J] = max(dataout,[],2);
    [~, J1] = max(target,[],2);
    counter = counter + length(find(J==J1));
    
    % expected log-likelihood error?
    err_cr = err_cr - sum(sum(target.*log(dataout)));
    %%%% wtf is there a minus instead of a plus??  It's this way in GEH's
    %%%% original code.
end

train_err = (Ncases*Nbatches - counter);
train_crerr = err_cr/Nbatches;

end