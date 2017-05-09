function [visDstrbs,visNums,inputDstrbs,inputNums,inputInd,...
    vishid,vishidinc,visbiases,visbiasinc,algorithm,Ncdsteps,...
    mwfunc,mvbfunc,bwfunc,bvbfunc,kwfunc,kvbfunc,...
    mw,mvb,bw,bvb,kw,kvb] = recurrentize(visDstrbs,visNums,...
    hidDstrbs,hidNums,vishid,vishidinc,visbiases,visbiasinc,epoch,params)
% recurrentize  Alter EFH parameters to train as a recurrent network


%-------------------------------------------------------------------------%
% Created: 12/26/16
%   by JGM
%-------------------------------------------------------------------------%


% set the recurrent algorithm with the requested number of CD steps
algorithm = params.algorithm;
Ncdsteps = params.Ncdsteps;

% for the (R)TRBM, store the old visible-unit properties
inputDstrbs = visDstrbs;
inputNums = visNums;
inputInd = sum(hidNums) + 1;

if ~strcmp(algorithm,'EFH')
    
    % augment the list of units in the visible layer
    visDstrbs   = [hidDstrbs, visDstrbs];
    %%%visDstrbs   = [visDstrbs, visDstrbs];
    %%%
    visNums     = [hidNums,   visNums];
    
    % augment the weight matrix, visible biases, and their increments
    %%% might be nice to synchronize this with reinitializeEFH.m....
    vishid = cat(1,0.005*randn([sum(hidNums),sum(hidNums)],...
        'like',vishid),vishid);
    vishidinc = cat(1,0*ones([sum(hidNums),sum(hidNums)],...
        'like',vishidinc),vishidinc);
    visbiases = cat(1,0*ones([sum(hidNums),1],...
        'like',visbiases),visbiases);
    visbiasinc = cat(1,0*ones([sum(hidNums),1],...
        'like',visbiasinc),visbiasinc);
    
    % augment the functions and parameters for the learning schedules
    mwfunc  = @(iEp)([params.mrw(iEp), params.mw(iEp)]);  mw  = mwfunc(epoch);
    mvbfunc = @(iEp)([params.mrb(iEp), params.mvb(iEp)]); mvb = mvbfunc(epoch);
    bwfunc  = @(iEp)([params.brw(iEp), params.bw(iEp)]);  bw = bwfunc(epoch);
    bvbfunc = @(iEp)([params.brb(iEp), params.bvb(iEp)]); bvb = bvbfunc(epoch);
    kwfunc  = @(iEp)([params.krw(iEp), params.kw(iEp)]);  kw = kwfunc(epoch);
    kvbfunc = @(iEp)([params.krb(iEp), params.kvb(iEp)]); kvb = kvbfunc(epoch);
else
    mwfunc  = @(iEp)(params.mw(iEp));       mw  = mwfunc(epoch);
    mvbfunc = @(iEp)(params.mvb(iEp));      mvb = mvbfunc(epoch);
    bwfunc  = @(iEp)(params.bw(iEp));       bw = bwfunc(epoch);
    bvbfunc = @(iEp)(params.bvb(iEp));      bvb = bvbfunc(epoch);
    kwfunc  = @(iEp)(params.kw(iEp));       kw = kwfunc(epoch);
    kvbfunc = @(iEp)(params.kvb(iEp));      kvb = kvbfunc(epoch);
end



end
