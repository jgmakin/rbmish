%% set up model: get samples of s and the means of the input layer
clear; clc;

Nsamples = 10000;
Ndims = 2;
Ngauss = 50;
Nbrnll = 20;

% "stimulus"
S = rand(Ndims,Nsamples);

% neuron means
M = rand(Ngauss,Ndims);
MEANgauss = M*S;

P(:,1) = rand(Nbrnll,1);
P(:,2) = 1 - P(:,1);
MEANbrnll = P*S;


% proj(Rgauss')*DECODEgauss = best guess for S
[U,D,Vgauss] = svd(MEANgauss');
lowDproj = U*D(:,1:2);
DECODEgauss = (lowDproj'*lowDproj)\lowDproj'*S';

% proj(Rbrnll')*DECODEbrnll = best guess for S
[U,D,Vbrnll] = svd(MEANbrnll');
lowDproj = U*D(:,1:2);
DECODEbrnll = (lowDproj'*lowDproj)\lowDproj'*S';


%% draw samples from the input layer

Ngauss = 2;
Rgauss = S;
% Rgauss = MEANgauss + randn(Ngauss,Nsamples);
Rbrnll = MEANbrnll > rand(size(MEANbrnll));

Rmat = cat(2,Rgauss',Rbrnll');
batchdata = shortdata(10,3,Rmat);



%% train an EFH
[numcases,numdims,numbatches] = size(batchdata);
params.t = Nbrnll;
params.N = numdims;
% params.numsUnits = [params.N params.N/2 params.N/2];
% params.typeUnits = {'GB','Bernoulli','Bernoulli'};
params.numsUnits = {params.N params.N*500};
params.typeUnits = {{'Bernoulli'},{'Bernoulli'}};
numsUnits = params.numsUnits;

% training params
params.NepochsMax = 90;                     % 50 for _Science_, but &c.
params.epsilonw =  200e-5;                  % 0.02 for all three
params.epsilonvb = 200e-5;
params.epsilonhb = 200e-5;
params.weightcost = 0.001;                  % 0.02 works well
params.initialmomentum = 0.5;
params.finalmomentum = 0.8;     % 0.9;      % 0.2 works well
params.counterMax = 8;
params.numtestbatches = 40;
params.Ncdsteps = 1;

% backprop params
params.BPmaxepoch = 5;
params.max_iter = 3;                        % number of linesearches
params.numMinibatches = 10;


% train
numRBMs = length(params.numsUnits)-1;
for iRBM = 1:numRBMs
    fprintf(1,'Pretraining Layer %i w/EFH: %d-%d \n',iRBM,numvis,numhid);
    restart = 1;
    
    % train
    tic; EFH; toc;
    
    % pack together weights for saving (hid => recog., vis => gener.)
    wts2{iRBM} = [vishid; hidbiases];
    wts2{numRBMs*2-iRBM+1} = [vishid'; visbiases];
    
    % for next time through
    batchdata = batchposhidmeans;   
    
end

    
%% compute the errors

% now test the model
Dout = updownDBN(Rmat,wts2,params,'means','quiet');


ShatgaussIN = Rgauss';
% projRgauss = Rgauss'*Vgauss;
% ShatgaussIN = projRgauss(:,1:2)*DECODEgauss;
figure(135);
subplot(1,2,1);
scatter(S(1,:),ShatgaussIN(:,1));
% axis equal
subplot(1,2,2);
scatter(S(2,:),ShatgaussIN(:,2));
% axis equal
EgaussIN = S - ShatgaussIN';


ShatgaussOUT = Dout(:,1:Ngauss);
% projRgauss = Dout(:,1:Ngauss)*Vgauss;
% ShatgaussOUT = projRgauss(:,1:2)*DECODEgauss;
figure(136);
subplot(1,2,1);
scatter(S(1,:),ShatgaussOUT(:,1));
axis equal
subplot(1,2,2);
scatter(S(2,:),ShatgaussOUT(:,2));
axis equal
EgaussOUT = S - ShatgaussOUT';


MSEgaussIN = EgaussIN*EgaussIN'/(Nsamples-1);
MSEgaussOUT = EgaussOUT*EgaussOUT'/(Nsamples-1);


projRbrnll = Rbrnll'*Vbrnll;
ShatbrnllIN = projRbrnll(:,1:2)*DECODEbrnll;
figure(137);
subplot(1,2,1);
scatter(S(1,:),ShatbrnllIN(:,1));
axis equal
subplot(1,2,2);
scatter(S(2,:),ShatbrnllIN(:,2));
axis equal
EbrnllIN = S - ShatbrnllIN';


projRbrnll = Dout(:,(end-params.t+1):end)*Vbrnll;
ShatbrnllOUT = projRbrnll(:,1:2)*DECODEbrnll;
figure(138);
subplot(1,2,1);
scatter(S(1,:),ShatbrnllOUT(:,1));
axis equal
subplot(1,2,2);
scatter(S(2,:),ShatbrnllOUT(:,2));
axis equal
EbrnllOUT = S - ShatbrnllOUT';


MSEbrnllIN = EbrnllIN*EbrnllIN'/(Nsamples-1);
MSEbrnllOUT = EbrnllOUT*EbrnllOUT'/(Nsamples-1);

% print
[det(MSEgaussIN) det(MSEbrnllIN) det(MSEgaussOUT) det(MSEbrnllOUT)]
%  0.0089    0.0013    0.0004    0.0004
1 - [diag(MSEgaussIN)'./var(S'); diag(MSEbrnllIN)'./var(S');...
    diag(MSEgaussOUT)'./var(S'); diag(MSEbrnllOUT)'./var(S')]
%    -1.2779   -0.7317
%     0.6082    0.3486
%     0.6085    0.5327
%     0.6454    0.5580



%% plot the estimates
    

[Y I] = sort(S');

for j = 1:2
    
    figure(101);
    subplot(1,2,j)
    hold on;
    plot(ShatgaussIN(I(:,j),j) - Y(:,j));
    plot(ShatgaussOUT(I(:,j),j) - Y(:,j),'r');
    % plot(Y(:,j),'k');
    hold off;
    
    figure(102);
    subplot(1,2,j)
    hold on;
    plot(ShatbrnllIN(I(:,j),j) - Y(:,j));
    plot(ShatbrnllOUT(I(:,j),j) - Y(:,j),'r');
    % plot(Y(:,j),'k');
    hold off;

end



%%
close all;
clear wts;
% wts{1} = wts2{1};
% wts{2} = wts2{4};
 wts = wts2;

numcases = 40;
numbatches = Nsamples/numcases;
DgaussOUT = zeros(numcases,Ngauss,numbatches);
Din = shortdata(numcases,Rmat);
Sshort = shortdata(numcases,S');
for iBatch = 1:numbatches
    
    % Dgauss = Din(:,1:Ngauss,iBatch);
    Dgauss = repmat(mean(Rgauss,2)',numcases,1);
    Dbrnll = Din(:,(Ngauss+1):end,iBatch);
    DELTA = inf;
    
    while norm(DELTA,'fro') > 1e-7
        
           
        if 0
            figure(135);
            subplot(1,2,1);
            scatter(Sshort(:,1,iBatch),Dgauss(:,1));
            subplot(1,2,2);
            scatter(Sshort(:,2,iBatch),Dgauss(:,2));
            title(num2str(iBatch))
            pause(0.05)
        end
        
        
        if 0
            projRgauss = Dgauss*Vgauss;
            ShatgaussIN = projRgauss(:,1:2)*DECODEgauss;
            figure(135);
            subplot(1,2,1);
            scatter(Sshort(:,1,iBatch),ShatgaussIN(:,1));
            subplot(1,2,2);
            scatter(Sshort(:,2,iBatch),ShatgaussIN(:,2));
            title(num2str(iBatch))
            pause(0.05)
        end
        
        
        Dout = updownDBN([Dgauss Dbrnll],wts,params,'means','quiet');
        DELTA = Dout(:,1:Ngauss) - Dgauss;
        Dgauss = Dout(:,1:Ngauss);
        
        
        if 0
            figure(123);
            imagesc(Dgauss);
            title(num2str(iBatch));
        end

        
    end
    DgaussOUT(:,:,iBatch) = Dgauss;
    
end

%%

RgaussOUT = longdata(DgaussOUT)';

scatter(MEANgauss(1,:),RgaussOUT(1,:))
scatter(MEANgauss(2,:),RgaussOUT(2,:))


ShatgaussOUT = RgaussOUT';
% % projRgauss = RgaussOUT'*Vgauss;
% % ShatgaussOUT = projRgauss(:,1:2)*DECODEgauss;
figure(136);
subplot(1,2,1);
scatter(S(1,:),ShatgaussOUT(:,1));
axis equal
subplot(1,2,2);
scatter(S(2,:),ShatgaussOUT(:,2));
axis equal
EgaussOUT = S - ShatgaussOUT';
MSEgaussOUT = EgaussOUT*EgaussOUT'/(Nsamples-1);


1 - diag(MSEgaussOUT)'./var(S')
%%%  0.7116    0.5850
%%%  0.7162    0.6048
%%%  0.7308    0.6124
%%%  0.7425    0.6382

%%%  0.7141    0.7079

%% what made this work?  
%   (1) It helped to have S be drawn from two independent uniform dstrbs,
%   rather than a MVN dstrb (w/correlations).
%
%   (2) You *don't* need to scale the *marginal* variance of the Gaussian
%   samples.  What needs to be unity is the *conditional variance*,
%   p(r|v)---and therefore, *presumably*, p(r|s).  One almost wants to say,
%   then, that you should noise up the trajectory info, X, with
%   unity-variance noise, before training the EFH on it.  Hmm...
%
%   (3) You needed to draw the Bernoulli samples with noise, rather than
%   just using the mode, i.e., rounding.  This example intuitively brings
%   out why the mode is so bad: 
%
%       -Take the simple case where there are only 4 (2D) locations, but 
%       they are encoded many times in the Bernoulli layer; say, for
%       concreteness, that each of the four points provides the mean for
%       each of 10 Bernoulli units.  If these units were to be sampled
%       "noiselessly" (rounded), all 10 associated units would give the
%       same answer; whereas, under a sampling scheme, the proportion of
%       each one of the four possible pairs (b/c 2D: [0,0],[0,1],etc.) 
%       would together give a better sense of the underlying mean, and
%       therefore location s.
%
%       -Now suppose that no point is represented more than once; e.g., say
%       that there are 10,000 samples S, drawn from a pair of (independent)
%       uniform dstrbs.  There's no longer anything to average over, so
%       shouldn't it be a matter of indifference whether or not we use
%       samples or modes?  
%   
%       Well, suppose we assumed some spatial continuity, some local
%       smoothness.  Then averaging samples over a small "patch" in the
%       input space would outperform the use of modes. E.g., modes will
%       label the entire upper right corner of S space with [1,1]---all
%       those points will be "ambiguated."
%
%       But whence this spatial continuity?  From, I think, the fact that
%       the functions in the network are continuous.  This is a very
%       important point.
%
%   (4) You also fixed your file to use params.t properly, i.e. as the
%   number of Bernoulli units, rather than as the number of left-side
%   units.
%
%   (5) You messed around w/the number of units---did that matter?























