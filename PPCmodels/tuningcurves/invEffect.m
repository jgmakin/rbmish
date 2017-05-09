function invEffect
% invEffect
%   inverse effectiveness?


%-------------------------------------------------------------------------%
% Revised: 02/24/16
%   -to reflect new format for numsUnits, typeUnits
% Created: 07/24/12
%   by JGM
%-------------------------------------------------------------------------%


clear; clc; close all;

load StdAllGains.mat


% params
%%% Nbatches = 25;
hidDstrbs = params.typeUnits{2};
hidNums = params.numsUnits{2};
Nhid = sum(hidNums);
M = 40;




% create a matrix full of a single pos/th: at the center of joint space
oo = ones(1,params.Ndims);
thcntr = scalefxn(0.5*oo,0*oo,1*oo,params.roboparams.thmin,params.roboparams.thmax);
poscntr = FK2link(thcntr,params.roboparams,1);
S = [poscntr; thcntr]';

% generate a grid of all the gains
[S,Q] = getStimuliTiled(M^2,class(Vbar),params);
Q.G = getUniformGains(M,params);
R = params.getData(S,Q);

% linearpart = @(x,w,b)(x*w + repmat(b,size(x,1),1));
% Vbar = linearpart(R,wts{1}(1:end-1,:),wts{1}(end,:));
V = invParamMap(R,wts{1}(1:end-1,:),wts{1}(end,:),hidDstrbs,hidNums,params);



% look for super/sub//additivity over all neurons, all gains
n = 7;
offset = 600;
h1 = figure(48);
h2 = figure(43);
% h3 = figure(59);
for i = 1:size(V,2)
    
    
    % "superadditivity" => A(i,j) > B(i,j), i.e. C < 0.
    A = reshape(V(:,i),M,M);                            % \in [0,1]
    B = repmat(A(:,1),1,M) + repmat(A(1,:),M,1);        % \in [0,2]
    C = B - A;                                          % \in [-1,2]
    D = B./A;
    
    if i < (n^2)
        set(0,'CurrentFigure',h1)
        subplot(n,n,i);
        imagesc(A,[0,1]);
        axis off
        
        set(0,'CurrentFigure',h2)
        subplot(n,n,i);
        imagesc(C,[-1 2]);
        % imagesc((C<0),[0 1]);
        axis off;
        
%         set(0,'CurrentFigure',h3)
%         subplot(n,n,i);
%         imagesc(D);
%         axis off;
    end
    
    caca(i) = sum(sum(C < 0));
end


% is it the case that we get superadditivity "when the modality-specific
% influences are very weak"?  That means that ...




% % generate a uniform (in prop) grid of stimuli
S = getUniformStims(M,params);

% gains
params.swing = 0;
gains =  [1 0; 0 1; 1 1; 15 0; 0 15; 15 15];

% malloc
V = zeros(size(S,1),Nhid,size(gains,1));

% generate tuning curves at different gains
[pool,HADBEENCLOSED] = parallelInit;
%%% no parfor??
for i = 1:size(gains,1);
    %%%% broken %%%%
    R = newTunerData(S,gains(i,:),params);
    %%%%%%%%%%%%%%%%
    V(:,:,i) = invParamMap(R,wts{1}(1:end-1,:),wts{1}(end,:),...
        hidDstrbs,hidNums,params);
end
if HADBEENCLOSED, delete(pool), end

% plot!
plotTCs(S(:,:,2),cat(3,V(:,:,1:3),sum(V(:,:,1:2),3)),params);% small inputs
plotTCs(S(:,:,2),cat(3,V(:,:,4:6),sum(V(:,:,4:5),3)),params);% big inputs




end
%-------------------------------------------------------------------------%