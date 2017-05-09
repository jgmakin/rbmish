function dispErrCovs(MargErrorStats,M,params,varargin)
% DISPERRCOVS   Display marginal error statistics
%   DISPERRCOVS(ERRORSTATS,PARAMS,VARARGIN) plots the error covariances
%   contrained in ERRORSTATS at the 95% confidence interval, centered at
%   their respective means.  The VARARGIN is an index for the figure
%   windows so as to avoid overwriting existing figures, if necessary.

%-----------------------------------------------------------------------%
% Revised: 11/05/14
%   -even tho this function is supposedly deprecated: changed the call to
%   tikzBarGraph to accommodate its new form
% Revised: 03/11/14
%   -changed 1D version to print tags.names (they may not exist!)
% Revised: 02/07/13
%   -consolidated loops by expliciting turning off the legend for the
%   centers
% Revised: 09/13/11
%   -added help and revision info (such as it is)
% Created: ??/??/??
%   by JGM
%-----------------------------------------------------------------------%

Nmods = length(params.mods);
Ndims = params.Ndims;
Nstats = length(MargErrorStats);
if isempty(varargin), offset = 0; else offset = varargin{1}; end

if Ndims > 1
    
    for iMod = 1:Nmods
        
        str = params.mods{iMod};
        legendCell = {};
        % figure(offset+i); hold on;
        hh = figure(); hold on;
                     
        for iStat = 1:Nstats
            for jMod = 1:Nmods
                
                stats = MargErrorStats{iStat}{jMod};
                if (jMod==iMod) || strcmp(stats.tags.src,'single')
            
                    if any(isnan(stats.cov(:)))
                        fprintf('covariance has NaN entries --- jgm\n\n');
                    else
                        
                        % pick out the right color
                        clrs = getColor(stats.tags.name);
                        
                        % plot the covariance
                        h = error_ellipse(stats.cov,stats.mu,'conf',.95); % ,'style',':');
                        set(h,'LineWidth',1,'Color',clrs);
                        legendCell = getLegend(legendCell,stats.tags);
                        
                        % plot the standard error of the mean
                        h = error_ellipse(stats.cov/M,stats.mu,'conf',.95);
                        set(h,'LineWidth',1,'Color',clrs);
                        hAnnotation = get(h,'Annotation');
                        hLegendEntry = get(hAnnotation','LegendInformation');
                        set(hLegendEntry,'IconDisplayStyle','off')
                    end
                end
            end
        end
        axis equal;
        
        % give it a title and legend
        titlestr = ['Error Covariances, ',str,...
            ' Populations (u_max = ',num2str(params.gmax(iMod)),')'];
        %%% titlestr = ['Error Covariances, ',str,...
        %%%       ' Populations ($u_{\text{max}} = ',num2str(g),'$)'];
        %%% title(titlestr,'Interpreter','Latex','Fontsize',15);
        %%% legend(legendCell,'Interpreter','latex','Fontsize',15);
        title(titlestr,'Interpreter','none');
        labelAxesRBMishly(params.mods{iMod})
        legend(legendCell,'Interpreter','none');
        hold off
        
        % convert
        switch  params.machine
            case {'kobayashi-maru','CUPCAKE','Themistocles'}
                legend off, title('');
                matlab2tikzWrapper(['2DerrorStats',date],hh);
                close;
            otherwise
                fprintf('skipping tikzBarGraph on this machine -- jgm\n');
        end
        
    end
       
    
else
        
    
    % malloc
    mus = NaN(Nmods,Nstats);
    vars = NaN(Nmods,Nstats);
    clrs = NaN(Nmods,Nstats,3);
    
    % init
    % legendCell = {};
    xaxislabels = {};
    
    for iStat = 1:length(MargErrorStats)
        for iMod = 1:Nmods
            stats = MargErrorStats{iStat}{iMod};
            mus(iMod,iStat) = stats.mu;
            vars(iMod,iStat) = stats.cov;
            clrs(iMod,iStat,:) = getColor(stats.tags.name);
            clrNames{iMod,iStat} = getColorName(stats.tags.name);
            % legendCell = getLegend(legendCell,stats.tags);
        end
        xaxislabels = {xaxislabels{:}, stats.tags.name};
    end
    
    % mean square errors (what you'll save to TikZ)
    MSEs = vars + mus.^2;
    %%% I suppose you could hack up an R^2, pretending the emissions are
    %%% linear-Gaussian, etc. etc.
    
    % plot
    for iMod = 1:Nmods
        barstar(mus(iMod,:),'estimate',...
            'error bias',xaxislabels,[],...
            squeeze(clrs(iMod,:,:)),['Error biases, ',params.mods{iMod}],...
            [],[]);
        
        barstar(vars(iMod,:),'estimate',...
            'error variance',xaxislabels,[],...
            squeeze(clrs(iMod,:,:)),['Error variances, ',params.mods{iMod}],...
            [],[]);
        
        
        %%% add units to the ylabel: radians and ??
        % write this one to tikz file
        switch params.machine
            case {'kobayashi-maru','CUPCAKE','Themistocles'}
                goodInds = ~isnan(MSEs(iMod,:));
                MSEerrorBars = zeros(sum(goodInds),2);
                tikzBarGraph(1:sum(goodInds),MSEs(iMod,goodInds)',...
                    MSEerrorBars,0,xaxislabels(goodInds)',...
                    'estimator','mean square error','',...
                    clrNames(iMod,goodInds)',...
                    2,0.32,'overlapping',{},...
                    ['1DerrorStats',params.mods{iMod},date]);
            otherwise
                fprintf('skipping tikzBarGraph on this machine -- jgm\n');
        end
    end
        
    
    
end
% set(h,'LineWidth',2);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function legendCell = getLegend(legendCell,tags)

% init with boldface
% varstr = '\mathbf{ ';
varstr = '{';

% modality symbol
switch tags.mod
    case 'Hand-Position'
        varstr = cat(2,varstr,'\visls}');
    case 'Joint-Angle'
        varstr = cat(2,varstr,'\props}');
    case 'Gaze-Angle'
        varstr = cat(2,varstr,'\gazes}');
    case 'Efference-Copy'
        varstr = cat(2,varstr,'\ctrls}');
end

% superscript and hat
if strcmp(tags.src,'single')
    newstr = cat(2,'$\hat{',varstr,'}^0');
else % multiple
    newstr = cat(2,'$\hat{',varstr,'}^1');
end

% theoretical or empirical?
if isfield(tags,'name')
    nametag = ['_\text{',tags.name,'}'];
    newstr = cat(2,newstr,nametag);
else
    if strcmp(tags.epist,'theoretical')
        newstr = cat(2,newstr,'_{opt}');
    end
end

% estimator or error?
if strcmp(tags.var,'error')
    newstr = cat(2,newstr,'-',varstr);
end


newstr = cat(2,newstr,'$');

legendCell = {legendCell{:},newstr};

end
%-------------------------------------------------------------------------%



% switch str
%     case 'Hand-Position'
%         legend({'$\hat{{\bf \theta}}^0$','$\hat{{\bf x}}^0$',...
%             '$\hat{{\bf x}}^0_{opt}$','$\hat{{\bf x}}^1$',...
%             '$\hat{{\bf x}}^1_{opt}$'},'Interpreter','latex');
%     case 'Joint-Angle'
%         legend({'$\hat{{\bf x}}^0$','$\hat{{\bf \theta}}^0$',...
%             '$\hat{{\bf \theta}}^0_{opt}$','$\hat{{\bf \theta}}^1$',...
%             '$\hat{{\bf \theta}}^1_{opt}$'},'Interpreter','latex');
%     case 'Gaze-Angle'
%         legend({'$\hat{{\bf e}}^0$','$\hat{{\bf e}}^0_{opt}$',...
%             '$\hat{{\bf e}}^1$','$\hat{{\bf e}}^1_{opt}$'},...
%             'interpreter','latex');
%     otherwise
%         fprintf('bad label; check params.mod -- jgm\n');
% end



%%%%%%
% case cupcake...
% The dipshits at the Mathworks have broken their DVI outputter, and
% therefore your plots.  See:
%
% http://stackoverflow.com/questions/19435848/how-to-make-fonts-available-to-the-latex-interpreter-in-matlab-r2013a
%
%%%%%%




