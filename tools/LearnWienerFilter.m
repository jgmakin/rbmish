function [betas, Rsq] = LearnWienerFilter(Snoisy,Sclean,Nhistorybins,...
    lambda,fraction_train)
% LearnWienerFilter   Fit a Wiener Filter
%
% USAGE:
%   Rsq = LearnWienerFilter(Snoisy,Sclean,T,Nhistorybins,lambda);
%
% Given noisy (Snoisy, [...] 
%   

%-------------------------------------------------------------------------%
% Created: 08/27/2017
%   by JGM
%-------------------------------------------------------------------------%

% Ns
Nsamples_per_sequence = size(Sclean,3);

if fraction_train < 1
    Ntrain = floor(fraction_train*Nsamples_per_sequence);
    Snoisy_heldout = Snoisy(:,:,(Ntrain+1):end);
    Sclean_heldout = Sclean(:,:,(Ntrain+1):end);
    [Snoisy_heldout, Sclean_heldout] =...
        historicizeSampleData(Snoisy_heldout,Sclean_heldout,Nhistorybins);

    Snoisy = Snoisy(:,:,1:Ntrain);
    Sclean = Sclean(:,:,1:Ntrain);
end

% reshape inputs and outputs to accommodate history, then regress
[Snoisy, Sclean] = historicizeSampleData(Snoisy,Sclean,Nhistorybins);
betas = linregress(Snoisy,Sclean,'L2 regularize',lambda);

% if you want to report results on *held out* data
if fraction_train < 1
	Snoisy_res = Sclean_heldout - Snoisy_heldout*betas;
    SStot = sum((Sclean_heldout - mean(Sclean_heldout,1)).^2);
    SSerr = sum(Snoisy_res.^2,1);
    Rsq = 1 - SSerr./SStot;
end


end

