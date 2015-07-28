function [statsL,statsN] = unimodalStim(D,S,wts,stimstr,params)
% UNIMODALSTIM  Test unimodal input to a trained DBN
%  USAGE: 
%       unimodalStim(testbatchdata0,Stest,wts,params);
%       unimodalStim(batchdata0,Strain,wts,params);
%  UNIMODALSTIM takes a dataset D (either batchdata0 or testbatchdata0) and
%  the *corresponding* causal stimuli (xtrain or xtest, resp.) along with
%  the weights of a trained-up DBN and its params, and generates images and
%  statistics of the inputs and outputs.
%-------------------------------------------------------------------------%
% Revised: 08/03/10
%   -added params as an argument to updown
% Created: 07/13/10
%   by JGM
%-------------------------------------------------------------------------%

% params
if strcmp(stimstr,'prop');
    DEADINPUT = 'vis';
    % params.NS = 'vis';
elseif strcmp(stimstr,'vis')
    DEADINPUT = 'prop';
    % params.NS = 'prop';
else
    error('unrecognized stimulus modality! -- jgm');
end



% wipe out one input
if strcmp(DEADINPUT,'vis')
    halfdata = killinput(D,1,params);
elseif strcmp(DEADINPUT,'prop')
    halfdata = killinput(D,2,params);
else
    error('unrecognized input to kill--- jgm');
end

% up-down pass
[halfin,Si] = longdata(halfdata,S);
clear D halfdata x;  
[~, halfout] = updown(halfin,wts,params,'means');

% display the results
% griddisp(halfin,xi,wts,params,halfout);

% get statistics using normal (CoM) decoding
[statsL,statsN] = estStatsCorePP(Si,params,'CoM',halfin,halfout);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function halfdata = killinput(D,m,params)

halfdata = zeros(size(D));
for i=1:size(D,1)
    for j = 1:size(D,3)
        T = displayshape(D(i,:,j),params);
        % fooP = zeros(size(fooP));
        halfdata(i,:,j) = [T{1}(:)*(m-1); T{2}(:)*(2-m)]';
    end
end


end
%-------------------------------------------------------------------------%