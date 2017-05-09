function Rcrnts = RNNforwardpass(Inputs,HidInit,W,b,dstrbs,nums,params)
% RNNforwardpass    Forward pass of a recurrent neural networks
%
% USAGE:
%   Rcrnts = RNNforwardpass(Inputs,HidInit,W,b,FXN,params);
%
% Inputs is a (T x Nin) matrix of input vectors; HidInit is the (1 x Nhid)
% row vector of initial hiddens; W is the matrix that maps:
%
%   [hiddens(t-1),inputs(t)] -> hiddens(t);
%
% b is the bias vector for that affine transformation; dstrbs specifies the
% types of hidden unit (and therefore the appropriate element-wise
% nonlinearities), nums the numbers of each type; and params is the usual
% parameter structure.  The output Rcrnts has size (T x Nhids).
%
% NB!! This is hard-coded to have the "hidden" (recurrent) units on the
% left.

%-------------------------------------------------------------------------%
% Revised: 12/23/16
%   -changed so that the output is shifted forward one in time!!  That is,
%   the first entry is no longer HidInit, the first input to the RNN ("0"),
%   but the first *output* of the RNN ("1"); and the last entry is the last
%   output ("T").  This is what EFH wants anyway.
% Revised: 12/09/16
%   -precomputed the "input biases" (outside the loop through time)
%   -changed to accept only one trajectory at a time.  This means that this
%   function will basically *only* be used inside EFH.m--which is ok, since
%   EFHfilter has overlapping functionality.
% Revised: 02/23/15
%   -added argument numsUnits
% Revised: 12/01/15
%   -now expects Inputs and HidInit with a different orientation
% Created: 05/29/15 (happy b'day VMO)
%   by JGM
%-------------------------------------------------------------------------%

USEGPU = 0;

% Ns
Nhids = size(HidInit,2);

% precompute the "input biases," i.e. Wrr*input + b.
Wrr = W(1:Nhids,:);              % hiddens on left!
Wir = W((Nhids+1):end,:);
Binput = Inputs*Wir + b;

% is it faster to do this *off* the GPU?
if USEGPU
    Rcrnts = forwardpass(HidInit,Wrr,Binput,dstrbs,nums,params);
else
    Rcrnts = forwardpass(gather(HidInit),gather(Wrr),gather(Binput),...
        dstrbs,nums,params);
    if isa(Inputs,'gpuArray'), Rcrnts = gpuArray(Rcrnts); end
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function Rcrnts = forwardpass(RcrntsNow,Wrr,Binput,dstrbs,nums,params)

% malloc
[T,Nhids] = size(Binput);
Rcrnts = zeros(T,Nhids,'like',RcrntsNow);

for t = 1:T
    RcrntsNow = invParamMap(RcrntsNow,Wrr,Binput(t,:),dstrbs,nums,params);
    Rcrnts(t,:) = RcrntsNow;
end

end
%-------------------------------------------------------------------------%
