function V = tuner(wts,params)
% TUNER experimental determination of tuning curves
%   TUNER determines the tuning curves of the hidden layer of a DBN
%   described by the cell array of weights WTS by histogramming hidden
%   nodes' response to various stimuli.
%
%-------------------------------------------------------------------------%
% Revised: 06/21/14
%   -largely rewrote, mostly to accommodate the matrix version of GTrespfxn
%   (and therefore to get rid of parfor), but also to fix some bugs.  You
%   haven't fully tested this version!  But it does not appear to be in use
%   anyway.
% Revised: 08/31/11
%   -rationalized
% Revised: 11/29/10
%   -changed to accomodate x in *true* coordinates
% Created: 07/29/10
%   by JGM
%-------------------------------------------------------------------------%

% init
N = params.N;
Ndims = params.Ndims;
g = params.g;
patchmin = params.margin*ones(1,Ndims);
patchmax = patchmin + params.respLength;
inputUnitType = params.typeUnits{1};

nPoints = 20;
nSamples = 5000;
nOutputs = size(wts{length(wts)/2},2);

V = zeros(nOutputs,nPoints,nPoints);

% sample evenly across prop and unevenly across vis
th1 = linspace(params.thmin(1),params.thmax(1),nPoints);
th2 = linspace(params.thmin(2),params.thmax(2),nPoints);
for i = 1:nPoints
    tic
    for j = 1:nPoints      
        
        
        %%%%%%%% hard-coded for Nmods = 2, Ndims = 2 %%%%%%%%
        G = unifSmplAboutCenter([g,g],params.swing,nSamples);
        
        theseTh = repmat([th1(i); th2(j)],[1,nSamples]);
        thPatch = scalefxn(theseTh,params.thmin,params.thmax,...
            patchmin,patchmax)';
        D0(:,1:N^Ndims) = PPCencode(thPatch,G(:,1),inputUnitType,params);    
        
        theseX = repmat(FK2link([th1(i); th2(j)],params,1)',[1,nSamples]);
        xPatch = scalefxn(theseX,params.posmin,params.posmax,...
            patchmin,patchmax)';
        D0(:,(N^Ndims + 1):(2*N^Ndims)) =...
            PPCencode(xPatch,G(:,1),inputUnitType,params);    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % push through network
        probs = D0;
        for layer = 1:length(wts)/2
            HIDFXN = params.typeUnits{layer+1};
            probs = feedforward(probs,wts{layer}(1:end-1,:),...
                wts{layer}(end,:),HIDFXN,params);
        end
        % average and bin
        V(:,i,j) = mean(probs);
    end
    fprintf('i = %i;    ',i);
    toc
end



rfielddisp(V);

end