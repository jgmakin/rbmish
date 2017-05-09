function figmap = EFHdisp(figmap,iEFH,epoch,wts,params,varargin)
% (figmap,datagenargs,params,iEFH,varargin)
% EFHdisp   EFH display
%   EFHdisp displays--and calculates--useful information associated with
%   the currently training EFH.
%
%   USAGE:
%       figmap = EFHdisp(figmap,datagenargs,params);
%
%       figmap = EFHdisp(figmap,datagenargs,params,vishid,hidbiases,...
%           visbiases,vishidinc,pvisstates,phidstates,qvisstates,epoch);


%-------------------------------------------------------------------------%
% Revised: 05/31/16
%   -removed longdata from final case ("otherwise") for consistency with
%   other cases
% Revised: 11/30/15
%   -total rewrite, mostly to incorporate new, more useful plots
%   -brought testData in here as a persistent variable
% Revised: 06/03/15
%   -modernized, somewhat
% Revised: 09/18/14
%   -renamed rbmdisp.m -> EFHdisp.m
% Cribbed: 01/28/11
%   -from rbm.m
%   by JGM
%-------------------------------------------------------------------------%

persistent Rtest Stest Qtest

% initialize or plot?
if epoch==0
    if iEFH==1
        [figmap,Rtest,Stest,Qtest] =...
            populatePlotHandles(figmap,params);
    end
    % otherwise do nothing
else
    figmap = updateFigmap(figmap,iEFH,epoch,wts,params,...
        Rtest,Stest,Qtest,varargin{:});
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [figmap,Rtest,Stest,Qtest] = populatePlotHandles(figmap,params)
% NB!!  If params.swing is 100%, the decoding error computed on these data
% will start *increasing* after some point (e.g., epoch 20).  Don't be
% fooled!  It may well still be decreasing, which you can see by changing
% Rtest so that it is generated with (e.g.) params.swing = 0.  This has
% happened before.

basenum = 2010;

figKeys = keys(figmap);
for iKey = 1:length(figKeys)
    figStruct = figmap(figKeys{iKey});
    if any(figStruct.DISP)
        figStruct.fig = figure(basenum+iKey);
        title(figKeys(iKey));
        hold on;
        figStruct.plothandle = plot(NaN,NaN);
        hold off;
        figmap(figKeys{iKey}) = figStruct;
    end
end

if any(figmap('TestError').DISP)
    fprintf('\n\nGenerating TESTING data...\n');
    if checkGPUavailability
        dataclass = 'gpuArray'; 
    else
        dataclass = 'double';
    end
    [Rtest,Stest,Qtest] = params.getTestData(dataclass);
else
    Rtest = []; Stest = []; Qtest = [];
end

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function figmap = updateFigmap(figmap,iEFH,epoch,dbnwts,params,Rtest,...
    Stest,Qtest,varargin)

% init
setColors;
figKeys = keys(figmap);
fprintf('epoch: %4i: ',epoch);

for iKey = 1:length(figKeys)
    figStruct = figmap(figKeys{iKey});
    if figStruct.DISP(epoch)
        
        if usejava('desktop'), set(0,'CurrentFigure',figStruct.fig); end
        switch figKeys{iKey}
            case 'ReconError'
                
                if ~isfield(figStruct,'data'), figStruct.data = []; end
                figStruct.data(end+1) = gather(varargin{8});
                if usejava('desktop')
                    x = [1:(length(figStruct.DISP)*(iEFH-1)),...
                        (1:epoch) + length(figStruct.DISP)*(iEFH-1)];
                    tmp = [repmat(figStruct.DISP,[iEFH-1,1]);figStruct.DISP(1:epoch)];
                    x = x(tmp);
                    hold on;
                    set(figStruct.plothandle,'XData',x,'YData',...
                        figStruct.data,'color','m');
                    hold off;
                    %%%%%%
                    %axis([0,100,0,100]);
                    %%%%%%
                end
                fprintf('recon error = %6.4e; ',figStruct.data(end));
                figmap(figKeys{iKey}) = figStruct;
                %%%%%%%%%%%%
                %%% change recon error to also be at RBM 1????
                %%%%%%%%%%%%
                
            case 'TestError'
                vishid = varargin{1};
                hidbiases = varargin{2};
                visbiases = varargin{3};
                
                % gather wts %%%% hmmm....
                numRBMs = length(dbnwts)/2;
                dbnwts{iEFH} = [vishid; hidbiases];
                dbnwts{numRBMs*2-iEFH+1} = [vishid'; visbiases'];
                wts(1:iEFH) = dbnwts(1:iEFH);
                wts((iEFH+1):2*iEFH) = dbnwts(end-iEFH+1:end); %%%?
                
                if ~isfield(figStruct,'data'), figStruct.data = []; end
                figStruct.data(end+1) =...
                    gather(params.testEFH(Rtest,Stest,Qtest,wts));
                if usejava('desktop')
                    x = [1:(length(figStruct.DISP)*(iEFH-1)),...
                        (1:epoch) + length(figStruct.DISP)*(iEFH-1)];
                    tmp = [repmat(figStruct.DISP,[iEFH-1,1]);figStruct.DISP(1:epoch)];
                    x = x(tmp);
                    hold on;
                    set(figStruct.plothandle,'XData',x,'YData',...
                        figStruct.data,'color',EFHcolor);
                    hold off;
                end
                fprintf('test error = %6.3e; ',figStruct.data(end));
                figmap(figKeys{iKey}) = figStruct;
                
                testErrs = gather(figStruct.data);
                save([getdir('code'),'RBM/EFHtest.mat'],'testErrs');
                
            case 'Recons'
                if usejava('desktop')
                    iCase = 1; % just do the first of Ncases
                    if length(params.numsUnits{1}) > 1
                        visInds = (params.numsUnits{1}(1)+1):sum(params.numsUnits{1});
                    else
                        visInds = 1:params.numsUnits{1};
                    end
                    pvisstates = varargin{5};
                    qvisstates = varargin{7};
                    colormap(gray);
                    Ttrue = displayshape(pvisstates(iCase,visInds),params);
                    subplot(2,1,1);
                    PPCplot(cat(2,Ttrue{:}),params,'data');
                    Tconfab = displayshape(qvisstates(iCase,visInds),params);
                    subplot(2,1,2);
                    PPCplot(cat(2,Tconfab{:}),params,'confabulations');
                end
                
            case 'Hiddens'
                if usejava('desktop')
                    phidstates = varargin{6};
                    colormap(gray);
                    imagesc(phidstates);
                end
                
            case 'Weights'
                if usejava('desktop')
                    vishid = varargin{1};
                    
                    clf;
                    colormap(gray);
                    rows = 4; cols = 5; space = 0;
                    ax = getCustomAxesPos(rows,cols,space);
                    
                    if length(params.numsUnits{1}) > 1
                        visInds = (params.numsUnits{1}(1)+1):sum(params.numsUnits{1});
                    else
                        visInds = 1:params.numsUnits{1};
                    end
                    hidInds = ceil(size(vishid,2)*rand(rows*cols,1));
                    for iRow = 1:rows
                        for iCol = 1:cols
                            Twts = displayshape(vishid(visInds,...
                                hidInds((iRow-1)*rows+iCol)),params);
                            imagesc(cat(2,Twts{:}),'Parent',ax(iRow,iCol));
                            set(ax(iRow,iCol),'xtick',[],'ytick',[]);
                        end
                    end
                end
                
            case 'WeightNorm'
                if ~isfield(figStruct,'data'), figStruct.data = []; end
                figStruct.data(end+1) = gather(norm(varargin{1})); % vishid
                if usejava('desktop')    
                    x = [1:(length(figStruct.DISP)*(iEFH-1)),...
                        (1:epoch) + length(figStruct.DISP)*(iEFH-1)];
                    tmp = [repmat(figStruct.DISP,[iEFH-1,1]);figStruct.DISP(1:epoch)];
                    x = x(tmp);
                    hold on;
                    set(figStruct.plothandle,'XData',x,...
                        'YData',figStruct.data,'color','b');
                    hold off;
                end
                figmap(figKeys{iKey}) = figStruct;
                
            case 'WeightVelocityNorm'
                if ~isfield(figStruct,'data'), figStruct.data = []; end
                figStruct.data(end+1) = gather(norm(varargin{4})); % vishidinc
                if usejava('desktop')
                    x = [1:(length(figStruct.DISP)*(iEFH-1)),...
                        (1:epoch) + length(figStruct.DISP)*(iEFH-1)];
                    tmp = [repmat(figStruct.DISP,[iEFH-1,1]);figStruct.DISP(1:epoch)];
                    x = x(tmp);
                    hold on;
                    set(figStruct.plothandle,'XData',x,...
                        'YData',figStruct.data,'color','b');
                    hold off;
                end
                figmap(figKeys{iKey}) = figStruct;
                
            otherwise
                fprintf('unknown plotting fxn! -- jgm\n');
        end
    end
end
drawnow;
fprintf('\n');

end
%-------------------------------------------------------------------------%



% title('weights');

% if DISPMAP(4)
%     fprintf('max <W_i>: %f, min <W_i>: %f\n',...
%         max(mean(vishid,2)),min(mean(vishid,2)));
%     fprintf('max{hb}: %f, min{hb}: %f\n',...
%         max(hidbiases),min(hidbiases));
%     fprintf('max{vb}: %f, min{vb}: %f\n',...
%         max(visbiases),min(visbiases));
%     fprintf('max{d_c}: %f, min{d_c}: %f\n',...
%         max(max((negvisSfctStats))),min(min((qvisstates))));
%     fprintf('max{sum(d)}: %f, max{sum(d_c)}: %f\n',...
%         max(posvisact),max(negvisact));
%     fprintf('max{h*d}: %f\n',max(max(negprods)));
% end

