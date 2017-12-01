function [X, Y] = historicizeSampleData(X,Y,Nsamples_per_subseq)
% historicizeSampleData  Matrix of sequences to matrix of subseqences
%
% USAGE:
%   X = historicizeSampleData(X,[],Nsamples_per_subseq)
%
% Given a tensor X of size
%
%   (Nseq x Nfeatures x Nsamples_per_seq), 
%
% transform into a matrix of size 
%
%   (Nseq*(Nsamples_per_seq-Nsamples_per_subseq) 
%       x Nfeatures*Nsamples_per_subseq)
%
% by including all the subsequences of length Nsamples_per_subseq that
% overlap at all but one point.
%
% This is useful for, e.g., the Wiener filter, for which the output X is
% the input to a regression.  In that case one also wants to transform the
% output tensor Y, which this function will also do.
%
% NB that this function also pads X with a column of ones!!

%-------------------------------------------------------------------------%
% Created: 08/27/17
%   by JGM
%-------------------------------------------------------------------------%

% Ns
[Nseq, Nfeatures, Nsamples_per_seq] = size(X);

% construct indices for subsequences of X  (Nsamples_per_subseq x Nsubseq)
subseqinds = (1:Nsamples_per_subseq)' +...
    (0:(Nsamples_per_seq-Nsamples_per_subseq));

% extract subsequences, reshape -> 
%   (Nsubseq*Nseq x Nfeatures*Nsamples_per_subseq)
try 
    X = reshape(X(:,:,subseqinds), Nseq, Nfeatures, Nsamples_per_subseq, []);
    X = reshape(permute(X,[4,1,2,3]), [], Nfeatures*Nsamples_per_subseq);
catch
    X = gather(X);
    Y = gather(Y);
    X = reshape(X(:,:,subseqinds), Nseq, Nfeatures, Nsamples_per_subseq, []);
    X = reshape(permute(X,[4,1,2,3]), [], Nfeatures*Nsamples_per_subseq);
end
X = [X, ones(size(X,1),1)];
    

% also reshape the output accordingly
if ~isempty(Y)
    Y = reshape(permute(Y(:,:,Nsamples_per_subseq:end),[3,1,2]),...
        [], size(Y,2));
end
    

end















