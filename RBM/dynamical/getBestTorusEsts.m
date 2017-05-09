function [ShatBest,eBest] = getBestTorusEsts(Shat,Snotwrapped,params)
% getBestTorusEsts  Minimal estimates and errors for toroidal Shat
% 
% USAGE:
%
%   [ShatBest,eBest] = getBestTorusEsts(Shat,Snotwrapped,mods,params)
%
% It is assumed that:
%
%   X -> Snotwrapped -> Swrapped -> R -> Shat -> ShatBest
%
% where the first arrow was accomplished with latent2stimFxns, the second
% consists of wrapping onto an N-torus, the third to encoding, and the
% fourth to decoding.  This function implements the fifth arrow, i.e. gets
% the best (closest) estimate to, and the corresponding error for, the
% original, non-wrapped stimuli S.  Except for params, all inputs and
% outputs have size:
%
%   (Nexamples x Ndims x Nmods)
%
%
% How does the function work?  First, the errors between Shat and Swrapped 
% are calculated.  Then these errors are multiplied by computing all 
% possible shifts by the range of the torus.  Since all errors will be 
% within *two* ranges, this set contains the minimal error, which is simply
% extracted with min.  Finally, ShatBest is created by adding this minimal 
% error back onto S, the Snonwrapped.
%
% NB: This function should *only* be called when params.walls = 'wrapping'.

%-------------------------------------------------------------------------%
% Revised: 01/05/17
%   -changed to make compatible with static data--part of the Grand
%   Revision.
% Revised: 08/29/14
%   -now runs on both mods at once
%   -returns Y as a 3-tensor (Ncases x Ndims*Nmods x T)
% Revised: 07/07/14
%   -replaced loop calculation of CZ with tensor calculation
% Created: ??/??/14
%   by JGM
%-------------------------------------------------------------------------%


% transform latents into stimuli, possibly wrapping
[Swrapped,srange] = wrapStimuli(Snotwrapped,params.smin,params.smax,params.N);
[Nexamples,Ndims,Nmods] = size(Swrapped);

% compute minimal errors
eAct = Shat - Swrapped;                 % Nexamples x Ndims x Nmods
eBck = eAct - srange;                   % "
eFwd = eAct + srange;                   % "
e = cat(4,eBck,eAct,eFwd);              % Nexamples x Ndims x Nmods x 3
[~,minInds] = min(abs(e),[],4);         % 
[indDim1,indDim2,indDim3] = ndgrid(1:Nexamples,1:Ndims,1:Nmods);
eBest = e(sub2ind(size(e),indDim1,indDim2,indDim3,minInds));

% "correct" Shat
ShatBest = Snotwrapped + eBest;
    
end


