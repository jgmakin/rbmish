function plotFilterErr(t0,LDSdata,params,varargin)

%-------------------------------------------------------------------------%
% Revised: 01/27/15
%   -removed internal getColor, changed to use external getColor.m
% Revised: 12/16/13
%   -X -> S etc.
% Created: ??/??/??
%   -by JGM
%-------------------------------------------------------------------------%



% init
samples = t0:1000;
Ndims = params.Ndims;
Ncases = params.Ncases;
Nmods = params.Nmods;
Nstates = size(params.dynamics.A,2);
smin = params.smin;
smax = params.smax;
if Ndims==2, th0 = get2DOutline(params.thmin,params.thmax,42); end
% NSind = strcmp(params.mods,params.NS);
% [Ncases,q,Nbatches] = size(x0);
Nfltrs = length(varargin);
Z = LDSdata.Z;

        
% malloc
thisShat = NaN(Ndims,length(samples),Nfltrs);
thisCvrn = NaN(Ndims,Ndims,length(samples),Nfltrs);

% loop
for iCase = 1:Ncases
    for iMod = 1:Nmods

        switch iMod
            case 1
                thisS(1:Ndims,1:length(samples)) = params.dynamics.C*...
                    squeeze(Z(iCase,1:Nstates,samples));
            case 2
                thisS(1:Ndims,1:length(samples)) = params.dynamics.H*...
                    squeeze(Z(iCase,Nstates+1:end,samples));
        end
        % thisS(1:Ndims,1:length(samples)) = squeeze(S0(iCase,:,iMod,samples));
        
        inds = 1+(iMod-1)*Ndims:(iMod*Ndims);
        for iFltr = 1:Nfltrs
            thisShat(:,:,iFltr) =...
                squeeze(varargin{iFltr}.Xpct(iCase,inds,samples));
            thisCvrn(:,:,:,iFltr) =...
                squeeze(varargin{iFltr}.Cvrn(iCase,inds,inds,samples));
        end
        
        
        % just the theta1 errors
        figure(67+iMod); clf;
        for iDim = 1:Ndims
            
            
            % top subplots
            subplot(2,Ndims,iDim); hold on;
            for iFltr = 1:Nfltrs
                thisE = thisShat(iDim,:,iFltr) - thisS(iDim,:);
                clr = getColor(varargin{iFltr}.name);
                plot(samples,thisE,'Color',clr);
                lgndstr{iFltr} = [varargin{iFltr}.name,' (',...
                    num2str(mean(thisE.^2)*10^3,'%0.02f'),')'];
            end
            legend(lgndstr)
            tstr = ['$\hat\prop_',num2str(iDim),'-\prop_',num2str(iDim),'$ vs. time'];
            plot(samples,zeros(size(samples)),'c');
            hold off;
            
            
            % bottom subplots
            subplot(2,Ndims,iDim+Ndims); hold on
            for iFltr = 1:Nfltrs
                stddevP = squeeze(sqrt(thisCvrn(iDim,iDim,:,iFltr)));
                clr = getColor(varargin{iFltr}.name);
                shadedErrorBar(samples,thisShat(iDim,:,iFltr),stddevP,...
                    {'Color',clr},1);
            end
            plot(samples,thisS(iDim,:),'c');
            %%% legend('raw','RBM','opt','true')
            
            plot(samples,smin(iDim,iMod)*ones(size(samples)),'k');
            plot(samples,smax(iDim,iMod)*ones(size(samples)),'k');
            plot(samples,(smin(iDim,iMod)+0.1)*ones(size(samples)),'k:');
            plot(samples,(smax(iDim,iMod)-0.1)*ones(size(samples)),'k:');
            
            tstr = ['$\prop_',num2str(iDim),'$ vs. time'];
            hold off;
            
        end
        % title(num2str(j))
        
        if Ndims == 2
            figure(6+iMod); clf; hold on;
            for iFltr = 1:Nfltrs
                clr = getColor(varargin{iFltr}.name);
                plot(thisShat(1,:,iFltr),thisShat(2,:,iFltr),'Color',clr);
            end
            plot(thisS(1,:),thisS(2,:),'c')
            % axis([params.thmin(1) params.thmax(1) params.thmin(2) params.thmax(2)]);
            plot(th0(:,1),th0(:,2),'k')
            axis equal
            title(num2str(iCase))
            hold off
        end
        
    end
    pause;
        
end

end