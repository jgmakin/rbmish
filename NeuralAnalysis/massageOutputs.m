function [Y,highfiringneurons,mu] = massageOutputs(R)
% find high-firing, zscore, transpose---etc.?

%-------------------------------------------------------------------------%
% Cribbed: 05/16/13
%   -from KF4HHS
%   by JGM
%-------------------------------------------------------------------------%


%%% every one of these transformations, except zscoring, hurt the
%%% crossvalidated MSE

highfiringneurons = (mean(R)>0);            % greater than 0.5 Hz
% R2 = sqrt(R(:,highfiringneurons));          % "variance stabilizing xform"
% R2 = R(:,highfiringneurons);
% Y = zscore(R)';
mu = mean(R);
Y = (R - repmat(mu,size(R,1),1))';


end