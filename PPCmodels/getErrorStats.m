function eStats = getErrorStats(testPosteriors,S)
% getErrorStats     Get bias and error covariance
%
% USAGE:
% 	eStats = getErrorStats(testPosteriors,S)
%
% Given a structure testPosteriors with a field called Xpct of size
% (Nexamples x Ndims x Nposteriors), computer the error with the "stimulus"
% tensor S (Nexamples x Ndims), and store its (nan)mean and (nan)covariance
% in the fields .Xpct and .Cvrn, resp., of a structure called eStats.
%
% It also writes a field .tags into the structure eStats, which is itself a
% structure array.
% array.

%-------------------------------------------------------------------------%
% Revised: 07/15/14
%   -changed from structure array to structure; i.e., put the different
%   posteriors into the last dimensions of the fields .Xpct and .Cvrn,
%   rather than different elements of an array eStats.
%   -made eStats.tags into a structure array
%   -added field .N
% Cribbed: 07/15/14
%   -from mastertest.m (see for older version history)
%   by JGM
%-------------------------------------------------------------------------%

[~,Ndims,Nposteriors] = size(testPosteriors.Xpct);

% malloc
eStats.Xpct = NaN(Ndims,Nposteriors,'like',S);
eStats.Cvrn = NaN(Ndims,Ndims,Nposteriors,'like',S);

for iPosterior = 1:Nposteriors
    e = testPosteriors.Xpct(:,:,iPosterior) - S;
    eStats.Xpct(1:Ndims,iPosterior) = nanmean(e,1);
    eStats.Cvrn(1:Ndims,1:Ndims,iPosterior) = nancov(e,1);
    eStats.N(iPosterior) = sum(~isnan(e(:,1)));
    
    eStats.tags(iPosterior).name = testPosteriors.srcs{iPosterior};
    %%%% this isn't necessarily true.....
    if iPosterior < (Nposteriors-1), src='single'; else src='multiple'; end
    eStats.tags(iPosterior).src = src;
end


end
