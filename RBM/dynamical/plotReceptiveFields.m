function plotReceptiveFields(filterdata,V0)
% plot some receptive fields of the dynamical RBM


%-------------------------------------------------------------------------%
% Revised: 12/18/13
%   -put X into "canonical structure" (a dimension for Nmods)
% Revised: 11/04/13
%   -functionized
% Created: 11/01/13
%   by JGM
%-------------------------------------------------------------------------%


% collect into useful form
[Ncases,Nneurons,T] = size(V0);
X = [];
for t = 1:T
    X = cat(4,X,filterdata(t).states);
end


% plot1DReceptiveFields(X,V0)
plot2DReceptiveFields(X,V0,1,2);
plot2DReceptiveFields(X,V0,1,3);
plot2DReceptiveFields(X,V0,1,4);
plot2DReceptiveFields(X,V0,2,4);






end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function plot1DReceptiveFields(X,V0)

% params
[Ncases,Nneurons,T] = size(V0);
Nstates = size(X,2);
nBins = 10;


% loop through neurons and dimensions
for iNeuron = 1:Nneurons
    for iState = 1:Nstates
        
        % collect all of data for this state
        foo = squeeze(X(:,iState,1,:));
        stims = foo(:);
    
        % collect all of the "firing rates" of this neuron
        foo = squeeze(V0(:,iNeuron,:));
        rates = foo(:);     
        
        % discretize stimulus space
        edges = linspace(min(stims)-5*eps,max(stims)+5*eps,nBins+1);
        [N, binInds] = histc(stims,edges);
        N = N(1:end-1);
        locs = repmat(binInds,1,nBins)==repmat(1:nBins,length(binInds),1);
        
        % compute average "rates" on the discretized line
        AvgRate = locs'*rates./N;
        sqddev = (rates - AvgRate(binInds)).^2;
        VarRate = locs'*sqddev./(N-1);
        SEMRate = sqrt(locs'*sqddev./(N.*(N-1)));
        
        
        % plot
        title(num2str(iNeuron));
        subplot(2,2,iState)
        bincenters = edges(1:end-1) + diff(edges)/2;
        shadedErrorBar(bincenters,AvgRate,1.96*SEMRate,'-r',1);
        thisaxis = axis;
        axis([thisaxis(1),thisaxis(2),0,1]);
    end
    pause()

end


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function plot2DReceptiveFields(X,V0,dimA,dimB)
% dimA specifies the *row* of AvgRate, dimB specifies the *column*.  And
% yes, showrfs is correctly flipping (across a horizontal axis) the data,
% so that the lower left corner of the plot has the lowest bin (rather
% than the upper left corner having it).


% N's
[Ncases,Nneurons,T] = size(V0);
% nBins = 10;
nBins = 8;
Nmax = 100;
%%% [~,randinds] = sort(rand(Nneurons,1));
%%% inds = randinds(1:Nmax);
inds = 1:Nmax;


% discrecretize space of stimulus A
foo = squeeze(X(:,dimA,1,:));
stimsA = foo(:);
edgesA = linspace(min(stimsA)-5*eps,max(stimsA)+5*eps,nBins+1);
[~, binIndsA] = histc(stimsA,edgesA);
%%% bincentersA = edgesA(1:end-1) + diff(edgesA)/2;

% discrecretize space of stimulus B
foo = squeeze(X(:,dimB,1,:));
stimsB = foo(:);
edgesB = linspace(min(stimsB)-5*eps,max(stimsB)+5*eps,nBins+1);
[~, binIndsB] = histc(stimsB,edgesB);
%%% bincentersB = edgesB(1:end-1) + diff(edgesB)/2;


% malloc
A = zeros(nBins*nBins,Nmax);

% loop through all neurons
for k = 1:length(inds)
    
    % which neuron?
    iNeuron = inds(k);
    
    % collect all of the "firing rates" of this neuron
    foo = squeeze(V0(:,iNeuron,:));
    rates = foo(:);
    
    % "flatten" the two-component subscripts into a set of single indices
    linInds = sub2ind([nBins,nBins],binIndsA,binIndsB);
    locsAB = double(repmat(linInds,1,nBins^2)==repmat(1:nBins^2,length(linInds),1));
    NAB = sum(locsAB)';
    AvgRate = reshape(locsAB'*rates./NAB,nBins,nBins);
    
    %%% careful! some squares in your 2D grid are never visited, so you
    %%% sometimes divide by zero here
    %%% isnan(AvgRate)
    %%% A HACK
    % AvgRate(isnan(AvgRate) & (reshape(locsAB'*rates,10,10)==0))=0;
    %%%
    
    
    sqddev = (rates - AvgRate(linInds)).^2;
    VarRate = reshape(locsAB'*sqddev./(NAB-1),nBins,nBins);
    SEMRate = reshape(sqrt(locsAB'*sqddev./(NAB.*(NAB-1))),nBins,nBins);
    
    
    % store these average rates
    A(:,k) = AvgRate(:);
    
    
end

h = figure;
showrfs(A,'black')
title([num2str(dimA),' vs. ',num2str(dimB)])
%%% as in "x vs. t," or "y vs. x"


fileloc = 'C:\Users\makin\Documents\#texs\ucsf\figs\dynamical\';
figtitle = ['RFs',num2str(dimA),num2str(dimB)];
fileext = '90.eps';
filename = [fileloc,figtitle,fileext];
saveTightFigure(h,filename)

end
%-------------------------------------------------------------------------%















