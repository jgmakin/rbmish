function V = tuner(wts,params)
% TUNER experimental determination of tuning curves
%   TUNER determines the tuning curves of the hidden layer of a DBN
%   described by the cell array of weights wts by histogramming hidden
%   nodes' response to various stimuli.
%
%-------------------------------------------------------------------------%
% Revised: 08/24/16
%   -changed to accommodate new version of PPCencode (which now does the
%   rescaling of stimuli within it, rather than without it).  This function
%   still has not been tested recently.
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
patchmin = params.margin*ones(1,Ndims);
patchmax = patchmin + params.respLength;
inputDstrbs = params.typeUnits{1};
inputNums = params.numsUnits{1};

% Ns
Npoints = 20;
Nsamples = 5000;
Noutputs = size(wts{length(wts)/2},2);

% malloc
V = zeros(Noutputs,Npoints,Npoints);

% sample evenly across prop and unevenly across vis
th1 = linspace(params.roboparams.thmin(1),params.roboparams.thmax(1),Npoints);
th2 = linspace(params.roboparams.thmin(2),params.roboparams.thmax(2),Npoints);
for i = 1:Npoints
    tic
    for j = 1:Npoints      
        
        
        %%%%%%%% hard-coded for Nmods = 2, Ndims = 2 %%%%%%%%
        G = UniformNormalDiracSampler(zeros(size(params.gmin),yrclass),...
            Inf,Nsamples,params.gmin,params.gmax,0);
        
        theseTh = repmat([th1(i) th2(j)],[Nsamples,1]);
        thPatch = scalefxn(theseTh,params.roboparams.thmin,...
            params.roboparams.thmax,patchmin,patchmax);
        D0(:,1:N^Ndims) = PPCencode(thPatch,G(:,1),...
            params.roboparams.thmin,params.roboparams.thmax,inputDstrbs,...
            inputNums,params);    
        
        theseX = repmat(FK2link([th1(i) th2(j)],params.roboparams,1),[Nsamples,1]);
        xPatch = scalefxn(theseX,params.roboparams.posmin,...
            params.roboparams.posmax,patchmin,patchmax);
        D0(:,(N^Ndims + 1):(2*N^Ndims)) = PPCencode(xPatch,G(:,1),...
            params.roboparams.thmin,params.roboparams.thmax,inputDstrbs,...
            inputNums,params);    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % push through network
        stdParams = D0;
        for layer = 1:length(wts)/2
            hidDstrbs = params.typeUnits{layer+1};
            hidNums = params.numsUnits{layer+1};
            stdParams = invParamMap(stdParams,wts{layer}(1:end-1,:),...
                wts{layer}(end,:),hidDstrbs,hidNums,params);
        end
        % average and bin
        V(:,i,j) = mean(stdParams);
    end
    fprintf('i = %i;    ',i);
    toc
end



rfielddisp(V);

end