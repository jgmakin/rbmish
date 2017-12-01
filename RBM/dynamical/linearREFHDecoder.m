function [wts,XhatStatic,RsqStatic,XhatDynamic,RsqDynamic] =...
    linearREFHDecoder(Rtrain,Xtrain,Qtrain,Rtest,Xtest,Qtest,wts,params,SStot)

%-------------------------------------------------------------------------%
% Cribbed: 09/29/17
%   -from filtersForNeuralData.m (JGM)
% Created: ~04/xx/16
%   by JGM
%-------------------------------------------------------------------------%


% if necessary, train an rEFH
if isempty(wts)
    
    dataclass = class(Rtrain);
    numsUnits = params.numsUnits;
    datagenargs = {};
    if isfield(params,'EACHBATCHISATRAJ')
        EACHBATCHISATRAJ = params.EACHBATCHISATRAJ;
    else
        EACHBATCHISATRAJ = 0;
    end
    if isfield(params,'Npretrain')
        Npretrain = params.Npretrain;
    else
        Npretrain = 0;
    end
    if isfield(params,'sparse')
        sparsitycost= params.sparse.cost;
        estimRate   = params.sparse.phidNewFrac;
        phidTarget  = params.sparse.phidTarget;
    else
        sparsitycost= 0;
    end
    getLatents  = params.getLatents;
    getData     = params.getData;
    NEFHs       = length(numsUnits)-1;
    Ncases      = params.Ncases;
    Nbatches    = params.Nbatches;
    Ntest       = params.NepochsMax + 1;
    NepochsMax  = params.NepochsMax;
    
    % init
    paramDisplay(params);
    wts = cell(NEFHs*2,1);
    
    % pretraining
    for iEFH = 1:NEFHs
        fprintf(1,'Pretraining Layer %i w/EFH: %d-%d \n',...
            iEFH,sum(numsUnits{iEFH}),sum(numsUnits{iEFH+1}));
        RESTART = 1;
        
        % train
        tic; EFH; toc;
        
        % pack together weights for saving (hid => recog., vis => gener.)
        wts{iEFH} = [vishid; hidbiases];
        wts{NEFHs*2-iEFH+1} = [vishid'; visbiases'];
    end
end

% test
[~,XhatStatic,RsqStatic,XhatDynamic,RsqDynamic] = testEFHBMI(...
    Rtest,Xtest,Qtest,wts,params,Rtrain,Xtrain,Qtrain,SStot);

end
