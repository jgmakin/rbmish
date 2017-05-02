function [R,Q] = getDataErlangMixture(S,Q,params)
% getDataErlangMixture
%
% USAGE:
%   R = getDataErlangMixture(S,Q,params)
% 
% For use with the model 'ErlangMixtureToy'.

%-------------------------------------------------------------------------%
% Cribbed: 01/02/17
%   from generateData.m
%   by JGM
%-------------------------------------------------------------------------%


% params
ks = params.shapeparams;
ths = params.scaleparams;
Nvis = sum(params.numsUnits{1});


S((sum(S,2)==0),end+1) = 1;
Th = repmat(shiftdim((S*[ks,ths])',-1),[Nvis/2,1]);
%%%%
%[~,randinds] = sort(rand(size(S,1)*Nvis/2,1));
%Th = reshape(shiftdim(Th,2),[size(Th,1)*size(Th,3),size(Th,2)]);
%Th = Th(randinds,:);
%Th = shiftdim(reshape(Th,[size(S,1),Nvis/2,2]),1);
%%%%
Th = reshape(Th,[Nvis,size(S,1)])';
%%% harmfunc = @(zz)(sum(1./(1:zz)));
%%% shapeparams = repmat(S*params.shapeparams,[1,Nvis]);
%%% scaleparams = repmat(S*params.scaleparams,[1,Nvis]);
%%% Mu1 = arrayfun(harmfunc,shapeparams-1) - double(vpa(eulergamma)) + log(scaleparams);
%%% Mu2 = shapeparams.*scaleparams;
%%% R = sampleT([Mu1,Mu2],'Erlang',params);
R = sampleT(Th,{'Erlang'},size(Th,2),params);



end