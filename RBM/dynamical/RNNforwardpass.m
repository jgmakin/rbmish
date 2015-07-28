function Hiddens = RNNforwardpass(Inputs,HidInit,W,b,FXN,params)
% RNNforwardpass    Forward pass of a recurrent neural networks
%
% USAGE:
%   Hiddens = RNNforwardpass(Inputs,HidInit,W,b,FXN,params);
%   
% Inputs is a (T x N) matrix of T input vectors; HidInit is the initial
% hidden-vector state; W is the matrix that maps:
%
%   [hiddens(t-1),inputs(t)] -> hiddens(t);
%
% b is the bias vector for that affine transformation; f is a string
% identifying the element-wise nonlinearity; and params is the usual params
% structure.
%
% NB!!! This is hard-coded to have the "hidden" (recurrent) units on the
% left (i.e., 'BP').

%-------------------------------------------------------------------------%
% Created: 05/29/15 (happy b'day VMO)
%   by JGM
%-------------------------------------------------------------------------%

% Ns
Nhids = size(HidInit,2);
[Ntraj,~,T] = size(Inputs);

% malloc/initialize
Hiddens = zeros(Ntraj,Nhids,T,'like',Inputs);
Hiddens(:,:,1) = HidInit;
for t = 1:(T-1)
    Augvectors = [Hiddens(:,:,t),Inputs(:,:,t)]; %%% hiddens on left!!
    Hiddens(:,:,t+1) = feedforward(Augvectors,W,b,FXN,params);
end



end