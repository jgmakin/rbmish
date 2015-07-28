function [f, df] = CG_func(wtsVec,Dim,varargin)
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
%-------------------------------------------------------------------------%
% Adapted: 6/02/10
%   -merged CG_Classif and CG_Encoder
% Adapted: 5/26/10
%   from the Salakhutdinov/Hinton code (CG_CLASSIFY.m)
%   by JGM
%-------------------------------------------------------------------------%

% check which objective fxn/model is being used
params = varargin{1};
data = varargin{2};
if strcmp(varargin{end},'CLASSIFY')
    targets = varargin{3};
    OBJ = 'CLASSIFY';
elseif strcmp(varargin{end},'ENCODE')
    targets = data;
    OBJ = 'ENCODE';
else
    error('unrecognized objective obj. for conjugate gradient -- jgm\n');
end

% params
numLayers = length(Dim) - 1;
N = size(data,1);

% unvectorize weights
offset = 0;
wts = cell(numLayers,1);
for i = 1:numLayers
    wts{i} = reshape(wtsVec(offset+1:offset+(Dim(i)+1)*Dim(i+1)),...
        Dim(i)+1,Dim(i+1));
    offset = offset + (Dim(i)+1)*Dim(i+1);
end

% propagate probabilities up, storing dE/dw's
[dataout dEdw tilde] = forwardprop(data,targets,wts,0,OBJ,params);

% (negative) expected log-likelihood (empirical expectation) for
if strcmp(OBJ,'ENCODE')
    % P(data|dataout) ~ Bernoulli (cross-entropy error)
    f = -sum(sum(targets.*log(dataout) + (1-targets).*log(1-dataout)))/N;
elseif strcmp(OBJ,'CLASSIFY')
    % multi-class data (cross-entropy error function)
    f = -sum(sum(targets.*log(dataout)))/N;
else
    error('unknown objective');
end

% vectorize dE/dw matrices (for minimize.m)
df = zeros(size(wtsVec));
offset = 0;
for i = 1:length(dEdw)
    len = length(dEdw{i}(:));
    df(offset+1:offset+len) = dEdw{i}(:);
    offset = offset + len;
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [output dEdw delta] = forwardprop(input,data,wts,layer,OBJ,params)
% FORWARDPROP   recursive backprop function
%
% Given the initial input (INPUT, DATA), a cell array WTS of the weight
% matrices (including the biases weights as the final column), and the
% initial LAYER (= 0), FORWARDPROP produces the final layer's output and
% a cell array of matrices dE/dw---one cell for each layer's (incoming)
% weights.  The output delta, used in the recursion, can be ignored (it's
% the final backprop delta).
% 
% NB: In accordance with Hinton's DBN, this network computes the errors of
% the "middle" layer (i.e., the top of the DBN) according to a linear
% rather than sigmoidal function. 
%-------------------------------------------------------------------------%
% Adapted: 4/27/10
%   -from the identically named fxn in CG_Encoder.  This version is suited
%       to classification rather than dimensionality reduction (this
%       amounts to removing the LINEAR units and accounting for the EXPON
%       ones).
% Created: 4/27/10
%   by JGM
%-------------------------------------------------------------------------%

% update the layer
layer = layer + 1;

% check the params for what kinds of nodes this layer has
m = length(params.numsUnits);
HIDFXN = params.typeUnits(m-abs(m-layer-1));    % (really)

% compute the activation in this layer via sigmoid or linear fxn
activ = feedforward(input,wts{layer}(1:end-1,:),wts{layer}(end,:),...
    HIDFXN,params);

% recursive backprop
if (layer < length(wts))                        % NON-OUTPUT LAYERS
    % relay final output
    [output dEdw delta] = forwardprop(activ,data,wts,layer,OBJ,params);          
    
    % compute dE/dw for outgoing weights
    activ = [activ ones(size(activ,1),1)];      % augment for biases
    dEdw{layer+1} = activ'*delta;
    
    % update 
    if strcmp(FXN,'Gaussian')
        delta = (delta*wts{layer+1}');
    else
        dydx = activ.*(1-activ);
        delta = (delta*wts{layer+1}').*dydx;    % update delta
    end
    delta = delta(:,1:end-1);                   % remove last col (??)
    
else                                            % OUTPUT LAYER
    output = activ;    
    delta = (output - data)/size(data,1);       % *** see (1) ***
    dEdw = cell(size(wts));                     % have to assign something
end

if (layer == 1)                                 % final (input) dE/dw
    dEdw{layer} = [input ones(size(input,1),1)]'*delta; 
end

end
%-------------------------------------------------------------------------%
% *** (1) ***
% In Hinton's original code, the output layer's delta for the
% CLASSIFICATION network doesn't include the normalization term; i.e., it
% has simply delta = (output - data).  I've included the normalization term
% and then likewise introduced it into the objective function f above.
%
% Why leave it out?? (since after all it's in the ENCODER version).
% Possibly the division creates smaller numbers and introduces round-off
% error---or, possibly, something worse.  Keep this in mind....
%-------------------------------------------------------------------------%