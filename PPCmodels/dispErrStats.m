function dispErrStats(eStats,modality)
% dispErrStats      Display error statistics
%   
%   USAGE:
%       dispErrStats(eStats,mod)
%
% Given the structure array eStats with fields:
%
%   xpct
%   Cvrn
%   tags.name
%
% plots the results in a useful way: a bar graph of MSEs for 1D errors, and 
% error-covariance ellipses (95% CI) centered at the mean error for 2D
% data.  The function also requires NSmod, the neutral-space modality, so 
% that it can label the axes appropriately.
%
% NB: This is the replacement for dispErrCovs.m

%-------------------------------------------------------------------------%
% Created: 07/15/14
%   -changed eStats to be a structure rather than structure array
%   -changed eStats.tags to be a structure array
% Created: 07/04/14
%   by JGM
%-------------------------------------------------------------------------%


[~, machine] = system('hostname');
machine = strtrim(machine);
[Ndims,Nstats] = size(eStats.Xpct);
%%% tags = [eStats.tags{:}];

switch Ndims
    case 1
        
        xaxislabels = {eStats.tags(:).name};
        for iStat = 1:Nstats
            clrNames{iStat} = getColorName(eStats.tags(iStat).name); 
        end
        MSEs = eStats.Xpct(:).^2 + eStats.Cvrn(:);
        %%% I suppose you could hack up an R^2, pretending the emissions 
        %%% are linear-Gaussian, etc. etc.
        
        switch machine
            case {'kobayashi-maru','CUPCAKE','Themistocles','domestica'}
                MSEerrorBars = zeros(sum(~isnan(MSEs)),2);
                goodInds = ~isnan(MSEs);
                tikzBarGraph(1:sum(goodInds),MSEs(goodInds),...
                    MSEerrorBars,0,xaxislabels(goodInds)',...
                    'estimator','mean square error','',...
                    clrNames(goodInds)',...
                    2,0.32,'overlapping',{},...
                    ['1DerrorStats-',modality,'-',date]);
            otherwise
                fprintf('skipping tikzBarGraph on this machine -- jgm\n');
                figure;
                bar(MSEs(~isnan(MSEs)));
                ylabel('MSEs');
                set(gca,'XTickLabel',xaxislabels(~isnan(MSEs)));
        end
        
    case 2
        
        for iStat = 1:Nstats
            clrNames{iStat} = getColorName(eStats.tags(iStat).name); 
        end
        
        switch  machine
            case {'kobayashi-maru','CUPCAKE','Themistocles',...
                    'zamfir','domestica'}
                
                % write out to tikz file
                tikzErrorEllipse(eStats,0.95,clrNames,...
                    ['2DerrorStats-',modality,'-',date]);
                
            otherwise
                
                fprintf('skipping tikz plot on this machine -- jgm\n');
                figure; hold on;
                
                for iStat = 1:Nstats
                    
                    xpct = eStats.Xpct(:,iStat);
                    Cvrn = eStats.Cvrn(:,:,iStat);
                    M = eStats.N(iStat);
                    if any(isnan(Cvrn(:)))
                        fprintf('covariance has NaN entries --- jgm\n\n');
                    else
                        
                        % pick out the right color
                        clrs = getColor(eStats.tags(iStat).name);
                        
                        % plot the covariance
                        plothandle = error_ellipse(Cvrn,xpct,'conf',.95); % ,'style',':');
                        set(plothandle,'LineWidth',1,'Color',clrs);
                        
                        % plot the standard error of the mean
                        plothandle = error_ellipse(Cvrn/M,xpct,'conf',.95);
                        set(plothandle,'LineWidth',1,'Color',clrs);
                        hAnnotation = get(plothandle,'Annotation');
                        hLegendEntry = get(hAnnotation','LegendInformation');
                        set(hLegendEntry,'IconDisplayStyle','off')
                    end
                    
                end
                labelAxesRBMishly(modality)
                axis equal;
        end
        
        
        
    otherwise
        error('you never programmed a case for Ndims > 2 -- jgm');
end


end
