function Yhat = forwardGenerate(Z0,t0,DISPLAY,lefthandinit,wts,params)

%-------------------------------------------------------------------------%
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

% params
Ncdsteps = 25; % 50;
[Ntraj,~,T] = size(Z0);
Nlayers = length(params.numsUnits);
Mlayers = length(wts)/2+1;
visDstrbs = params.typeUnits{Mlayers-1};
visNums = params.numsUnits{Mlayers-1};
hidDstrbs = params.typeUnits{Mlayers};
hidNums = params.numsUnits{Mlayers};

% useful weights
if Mlayers==Nlayers
    input0 = sum(hidNums) + 1;
    inputDstrbs = [hidDstrbs,visDstrbs];
    inputNums = [hidNums,visNums];
else
    input0 = 1;
    inputDstrbs = visDstrbs;
    inputNums = visNums;
end
Wuz = wts{Mlayers-1}(1:(input0-1),:);
Wyz = wts{Mlayers-1}(input0:end-1,:);
Wzv = wts{Mlayers}(1:end-1,:);
by  = wts{Mlayers}(end,input0:end)';
bv  = wts{Mlayers}(end,:)';

% malloc/init
thisU = Z0(:,1:(input0-1),t0); % i.e., all or none
Yhat = zeros(Ntraj,sum(params.numsUnits{1}),T,'like',Z0);

fprintf('Forward Gibbs sampling...')
tic
for t = (t0+1):T
    fprintf('.')
    
    % initialize right-half Gibbs sampling
    Yrand = rand(Ntraj,params.numsUnits{Mlayers-1});
    if DISPLAY
        figure(5413); clf;
        topo = displayshape(Yrand(2,:),params);
        imagesc(topo{1})
        pause;
    end
    
    % first (one-half-Gibbs-step) ultimate layer
    V = cat(2,thisU,Yrand);
    phidparams = invParamMap(V,wts{Mlayers-1}(1:end-1,:),...
        wts{Mlayers-1}(end,:),hidDstrbs,hidNums,params);
    phidstates = sampleT(phidparams,hidDstrbs,hidNums,params);
    %%% might want to use means!
    
    % clamp left half, Gibbs sample right half
    bz = wts{Mlayers-1}(end,:) + thisU*Wuz;
    [~, qhidstates] = CDstepper(phidstates,Wyz,by,bz,...
        hidDstrbs,visDstrbs,hidNums,visNums,Ncdsteps,params);
    
    % one more, with means
    Ybar = invParamMap(qhidstates,Wyz',by',visDstrbs,visNums,params);    
    Zbar = invParamMap(Ybar,Wyz,bz,hidDstrbs,hidNums,params);
    Vbar = invParamMap(Zbar,Wzv,bv',inputDstrbs,inputNums,params);
    Yf = Vbar(:,input0:end);
    
    % what will be the left-hand units at the next time step?
    switch lefthandinit
        case 'confabhidmeans'
            thisU = Vbar(:,1:(input0-1));
        case 'confabhidstates'
            %%% ? this will probably fail for DBNs....
            thisU = sampleT(Vbar(:,1:(input0-1)),hidDstrbs,hidNums,params);
            %%%
        case 'inferredhidmeans'
            thisU = Z0(:,1:(input0-1),t); % all or none
    end
    
    % now push down to the first layer
    for iLayer = (Mlayers+1):length(wts)
        jLayer  = Nlayers-abs(Nlayers-iLayer-1); % really!
        dstrbs  = params.typeUnits{jLayer};
        nums    = params.numsUnits{jLayer};
        Yf = invParamMap(Yf,wts{iLayer}(1:end-1,:),...
            wts{iLayer}(end,:),dstrbs,nums,params);
    end
    Yhat(:,:,t) = Yf;
    
end
toc
fprintf('\n')


end
