function Pi = getCatProbs(Th,DSTRB)
% getEFHcatProbs    Convert parameters to category probabilities
%
% USAGE:
%
%   Pi = getCatProbs(Th,DSTRB)
%
% For parameters Th corresponding to the distribution specified by DSTRB,
% convert to category probabilities.  For 'Categorical', Th contains the
% probabilities of the first N-1 categories, so the final-category prob. is
% simply concatentated.  For 'Bernoulli', a moore complicated conversion is
% implemented.


%-------------------------------------------------------------------------%
% Created: 02/19/16
%   by JGM
%-------------------------------------------------------------------------%

switch DSTRB
    case 'Bernoulli'
        if 2^size(Th,2) > 40
            fprintf('warning: you haven''t programmed in this case -- jgm');
            Pi = Th;
        else
        A = permute(decimalToBinaryVector((1:2^size(Th,2))-1),[3,2,1]);
        Pi = permute(prod(Th.*A + (1-Th).*(1-A),2),[1,3,2]);
        end
    case 'Categorical'
        Pi = [Th, 1-sum(Th,2)];
    otherwise
        Pi = Th;
        %%%error('you haven''t programmed in this case -- jgm');
        fprintf('warning: you haven''t programmed in this case -- jgm');
end

