function griddisp(D,x,wts,params,varargin)
% USAGE: griddisp(batchdata0,x0,wts,params);
%        griddisp(testbatchdata0,x0,wts,params);
% displays inputs and outputs of "autoencoder," vis and prop, in one plot

%-------------------------------------------------------------------------%
% Revised: 12/10/12
%   -changed to accommodate new form of shatL and shatN: matrices 
%   (Ndims x Nmods) rather than vectors.
% Revised: 08/20/12
%   -finally fixed the third subplot for Nmods=3; made associated changes
% Revised: 05/09/11
%  -eliminated subroutine crossmark, b/c all it did was call mark.m and
%  decode, but the decoding can now be retrieved with the estError call.
%   -and made some other little changes, e.g. incorporating params.NS
% Revised: 01/31/11
%   -changed the way estGather works, so crossmark.m got changed
% Revised: 01/11/11
%   -changed to accomodate griddisp1D
% Revised: 11/29/10
%   -changed makecross to use world2grid.  NB that this subtly changes the
%   results!!
% Revised: 07/28/10
%   -put updown into this fxn (so you don't have to run that first)
% Created: 06/24/10
%   by JGM
%-------------------------------------------------------------------------%

% init
Ndims = params.Ndims;
Nmods = params.Nmods;
smin = params.smin;
smax = params.smax;
k = find(strcmp(params.mods,params.NS));


% these are hard-coded for the usual params.NS
switch params.MODEL
    case {'2Dinteg','1Dinteg','MCDexperiment'}
        ttlstr = '(input,output): {\color{red}F^{-1}(vis)}, prop (bump)';
        cbnmat = [0 1; 1 0];
    case {'1Daddition', '2Daddition'}
        ttlstr = '(input,output): {\color{red}F(prop)}, vis+eye (bump)';
        cbnmat = [0 1 1; 1 0 -1; 1 -1 0];
    otherwise
        error('weird number of input modalities -- jgm');
end
fxn = @(shat)(shat*cbnmat(:,k));


% put in long format
[Di,xi] = longdata(D,x);
Nexamples = size(Di,1);

% check if the second data set has been provided via a varargin
if isempty(varargin)
    [~, Do] = updown(Di,wts,params,'means');
else
    Do = varargin{1};
end



% malloc
errL = zeros(4,1);              errN = zeros(4,1);
eL = zeros(Ndims,2*Nmods);      eN = zeros(Ndims,2*Nmods);

% loop
figure('Position',[100,50,1300,750]); colormap(gray);
for iExample = 1:Nexamples
    
    % plot inputs and outputs (via estError)
    % the estimates are marked in red, the true stimuli in green
    subplot(4,1,1);                             % inputs
    [eL(:,1:Nmods), eN(:,1:Nmods), ~, shatNi] =...
        estError(Di(iExample,:),squeeze(xi(iExample,:,:)),params,'CoM',1);
    subplot(4,1,2);                             % outputs
    [eL(:,Nmods+1:2*Nmods), eN(:,Nmods+1:2*Nmods), ~, shatNo] =...
        estError(Do(iExample,:),squeeze(xi(iExample,:,:)),params,'CoM',1);
    
    
    % plot the population code of the neutral-space variable
    subplot(4,1,3);
    Ti = displayshape(Di(iExample,:),params);
    To = displayshape(Do(iExample,:),params);
    PPCplot(cat(2,[Ti{k},To{k}]),params,ttlstr);
    
    % mark the spot on this plot given by decoding the other population(s)
    hold on
    shatIN = fxn(shatNi);
    shatOUT = fxn(shatNo);
    mark(shatIN,[smin(:,k),smax(:,k)],max(Ti{k}),0,'r',params);
    mark(shatOUT,[smin(:,k),smax(:,k)],max(To{k}),1,'r',params);
    hold off
    
    
    % make a bar plot of the errors
    subplot(4,1,4);
    for j = 1:4
        errL(j) = norm(eL(:,j));
        errN(j) = norm(eN(:,j));
    end
    % bar([errL(1) errL(3) errL(2) errL(4)]);
    bar([errL(2) errL(4)]);
    % set(gca,'XTickLabel',{'errV_in','errV_out','errP_in','errP_out'});
    set(gca,'XTickLabel',{'errP_in','errP_out'});
    
    % break?
    figure(1)
    reply = input('quit? y/n  ','s');
    if strcmp(reply,'y')
        close;
        fprintf('bye\n');
        break
    end
    
end

% colormap gray;


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


% vmin = params.posmin;
% vmax = params.posmax;
% patchmin = [params.margin; params.margin];
% patchmax = patchmin + [params.respLength; params.respLength];
% granularity = params.granularity;
% vPatch = scalefxn(v,vmin,vmax,patchmin,patchmax);    % pos -> stddev un.
% vGrid = round(pos*granularity);                   % stddev units -> grd
% units