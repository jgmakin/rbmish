function Xhat = catAssign(X,Xhat)
% catAssign     Assign categorized data to best true categories
% 
% USAGE:
%
%   Xhat = catAssign(X,Xhat)
% 
% For each of Nexamples, the user supplies:
%
%   Xhat (Nexamples x N): the prob. of each of N fictive categories;
%   X (Nexamples x M): the true assignment to one of M true categories.
%
% This function assigns model (or "fictive") categoties to data (or "true")
% categories, based on their correlations, with the additional requirement
% that all data categories be assigned *exactly one* fictive category.
% 
% It is impossible for a data category to correlate with no model category:
% for every trial (and therefore every data category), *some* model
% categorization will be made.  Thus, a function from data to model
% categories is guaranteed to exist.  Nevertheless, the inverse mapping may
% not exist (model to data), since the (forward) map can fail to be either
% injective or surjective.  Additional rules are therefore needed to ensure
% that the inverse map exists.
% 
% If (e.g.) two different data categories correlate most highly with the
% same model category, the data category with the higher correlation is
% assigned that model category.  The "losing" data category is then
% assigned the model category with which it correlates *second* most.  The
% procedure repeats recursively in the obvious way.  This prevents failures
% of injectivity.
% 
% If a model category is not ultimately assigned to a data category--and
% this is likely, since usually N > M--then that model category is simply
% discarded.  This is true even if that model category is, for some trials,
% the *most* likely.

%-------------------------------------------------------------------------%
% Created: 02/12/16
%   by JGM
%-------------------------------------------------------------------------%

C = Xhat'*X;                        % compute the inner product ("corr.")
[aa,bb] = max(C,[],2);              % find best corr. for each fictive cat
[~,cc] = sort(aa,'descend');        % get sorted fictive-cat inds
dd = bb(cc);                        % corr. inds sorted by corr. size
ee = dd(1:size(C,2));               % first M "
ff = cc(1:size(C,2));               % first M sorted fictive-cat inds
gg = ff(ee);                        % first M sorted fictice-cat inds, 
                                    %  sorted by first M corr. inds
Xhat = Xhat(:,gg);                  % select and assign categories

end