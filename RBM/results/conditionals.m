% conditionals
%   Computes the conditional errors for a few different directions.  This
%   file is called by standardCond, nonflatpriorCond.

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -cosmetic changes
% Revised: 05/24/11
%   -pulled out from standardCond.m
% Created: 05/20/11
%   by JGM
%-------------------------------------------------------------------------%

%%
% init
load(wtsfile);
% load(NNfile);
nDirections = 6;
params.swing = 0;
Ndims = params.Ndims;


cntr = scalefxn(0.5*ones(Ndims,1),zeros(Ndims,1),ones(Ndims,1),...
    params.thmin,params.thmax);
% inc = 8/150;
inc = 8/150;
% inc = 0.03;
% inc = 2/150;

% compute results
if RESTART
    ErrorStatsArray = [];
    for iDirection = 1:nDirections
        
        phi = 2*pi*iDirection/nDirections;
        R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
        % mu = R*[inc; 0] + 0.5;
        % [p0.mu p0.cov] = priorcreator(mu,0,params);
        p0.mu = R*[inc; 0] + cntr;
        p0.cov = zeros(Ndims);
        
        % ErrorStats = test(wts,params,'prior',p0);
        ErrorStats = test(wts,params,'prior',p0,...
            'propagation','Nsamples','numsamples',15);
%         [ErrorStats net] = nnDecode(wts,'prop',params,...
%             'pretrained',net,'prior',p0);
        ErrorStatsArray = [ErrorStatsArray; ErrorStats];
    end
    RESTART = 0;                    % just to be safe
else
    load(condfile);
end



%%
% plot
close all
for iDirection = 1:nDirections
    
    phi = 2*pi*iDirection/nDirections;
    R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
    mu = R*[inc; 0] + 0.5;
    [p0.mu,p0.cov] = priorcreator(mu,0,params);
    
    MargErrorStats = ErrorStatsArray(iDirection,:);
    for j = 1:length(MargErrorStats)
        for iMod = 1:params.Nmods
            MargErrorStats{j}{iMod}.mu = MargErrorStats{j}{iMod}.mu + p0.mu;
        end
    end
    
    dispErrCovs(MargErrorStats([1,3,4]),params.Ncases*Nbatches,params);
    
    
    deg = num2str(round(180*phi/pi));
    figname{iDirection} = ['cond' deg];
end



% %%
% xmax = 0;
% ymax = 0;
% for iDirection = 1:nDirections
%     h = figure(2*iDirection-1);
%     close(h)
%     
%     figure(2*iDirection)
%     axis tight
%     v = axis;
%     xx = v(2) - v(1);
%     yy = v(4) - v(3);
%     xmax = (xx > xmax)*xx + (xx <=xmax)*xmax;
%     ymax = (yy > ymax)*yy + (yy <=ymax)*ymax;
% 
% end
% 
% 
% for iDirection = 1:nDirections
%     
%     figure(2*iDirection)
%     v = axis;
%     v(2) = v(1) + xmax;
%     v(4) = v(3) + ymax;
%     axis(v);
%     if iDirection < nDirections
%         legend('off')
%     end
%     % saveas(gcf,['results\figs\',figname{iDirection}],'pdf');
% end
% 
