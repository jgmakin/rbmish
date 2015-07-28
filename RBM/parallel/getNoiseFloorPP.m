function stats = getNoiseFloorPP(D,x,wts,stimstr,params)
% [mu,covVo,covPo] = getNoiseFloorPP(D,x,wts,stimstr,params)
% D = (unmolested) input data
% x = corresponding (true) underlying stimuli
% wts = the RBM wts learned for these data
% stimstr = a string denoting which stimulus ('vis' or 'prop') to *keep*
% params = the global universal parameter structure
%
% TO DO: 
% (0) fxnize

clc; close all

% init
Nmods = params.Nmods;
Ndims = params.Ndims;
alpha = 0.001;
[Ncases,Nvis,Nbatches] = size(D);
% e_v = zeros(2,2,numBatches*numCases);
% e_p = zeros(2,2,numBatches*numCases);
% e = zeros(Ndims,2*Nmods,numBatches*numCases,Nmods);         % 2 b/c input *and* output
eN = zeros(Ncases*Nbatches,Ndims,Nmods);
eL = zeros(Ncases*Nbatches,Ndims,Nmods);

% zero out half the data
switch stimstr
    case 'prop'
        zvec = 1:Nvis/2;
    case 'vis'
        zvec = Nvis/2+1:Nvis;
        alpha = alpha/100000/100;
    otherwise
        error('unrecognized stimulus (should be ''vis'' or ''prop''');
end
D(:,zvec,:) = 0;

[Di,x] = longdata(D,x);
[avgErr,Do] = updownfast(Di,wts,params);

% loop
% matlabpool open local 8
% parfor (j = 1:size(Di,1))
for iExample = 1:size(Di,1)
    
    % init
    d_in = Di(iExample,:);
    d_out = Do(iExample,:);
    xtrue = x(iExample,:);

    %display?
    h = figure; colormap(gray);
    T = displayshape(d_in,params);
    PPCplot(cat(2,T{:}),params,'responses: r_{vis}, r_{prop}');
    saveas(h,['figs/xform',num2str(1),'.eps']);
    pause()
    
    % first iteration
    counter = 1;
    e0 = norm(d_in - d_out)/Nvis;
    e0_old = e0/(alpha-1)-1;
    
    % iterate with clamped half-input until convergence
    while ((e0_old - e0)/e0_old > alpha) && (counter < 100)
        
        % update
        d_in(zvec) = d_out(zvec);
        e0_old = e0;
        counter = counter + 1;
        
        % display?
        h = figure; colormap(gray);
        T = displayshape(d_out,params);
        PPCplot(cat(2,T{:}),params,'responses: r_{vis}, r_{prop}');
        drawnow
        saveas(h,['figs/xform',num2str(counter),'.eps']);
        pause()
        
        % compute error
        [avgErr,d_out] = updownfast(d_in,wts,params);
        e0 = norm(d_in - d_out)/Nvis;
    end
    if counter > 100
        fprintf('counter ran out!\n');
    end
    
    % compute error for both d_in and d_out
    [eL(iExample,:,:),eN(iExample,:,:)] = estError(d_out,xtrue,params,'CoM',0);
    
    % display stuff
    fprintf('counter exit: %i  example: %i\n',counter,iExample);
%     T = displayshape(d_out,params);
%     imagesc([T{1},T{2}]);
%     pause()
    % griddisp(d_out,wts,params,d_in);
    % pause()
end
% matlabpool close

% do it all on e_v?
statsL = cell(Nmods,n);
statsN = cell(Nmods,n);
for iMod = 1:Nmods
    statsL{iMod}.cov = nancov(eL(:,:,iMod));
    statsL{iMod}.mu = nanmean(eL(:,:,iMod));
    statsN{iMod}.cov = nancov(eN(:,:,iMod));
    statsN{iMod}.mu = nanmean(eN(:,:,iMod));
end



end