function [EMMparams,inds] = EFparams2ErlangParams(Pizq,Thyq)
% Standard (EF) params/per trial to Erlang-mixture-model parameters

%-------------------------------------------------------------------------%
% Cribbed: 06/02/16
%   from evaluateEMMEFH.m
%   by JGM
%-------------------------------------------------------------------------%

% find the *average* set of Erlang parameters associated w/each category
[~,bb] = max(Pizq,[],2);
inds = unique(bb);
MeanParamsByCat = arrayfun(@(iCat)(mean(Thyq(bb==iCat,:),1))',inds,...
    'UniformOutput',false);
% STDofParamsByCat = arrayfun(@(iCat)(sqrt(var(Thyq(bb==iCat,:))))',...
%     inds,'UniformOutput',false);
% CVofParamsByCat = [STDofParamsByCat{:}]./MeanParamsByCat;
% keyboard
MeanParamsByCat = [MeanParamsByCat{:}];
EMMparams.ks = MeanParamsByCat(1,:)';
EMMparams.mus = MeanParamsByCat(2,:)';
EMMparams.pis = mean(Pizq(:,inds))';

end