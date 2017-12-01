function filterRelativeContributions


%-------------------------------------------------------------------------%
% Created: 10/18/17
%   by JGM
%-------------------------------------------------------------------------%

monkey = 'Indy';

% parameters
binwidth = 64;
trainingtime = 320;
fracneuron = 1.0;
swept_param = 'binwidths';
obsv_dstrb = 'Poisson';

TOPLOT = 0;



for iSession = 21 %1:37
    % other useful things
    if checkGPUavailability, dataclass='gpuArray'; else, dataclass='double';end
    
    
    % load and test an rEFH
    [wts,params,model_full_path] = loadEFHspikecounts(monkey,iSession,...
        binwidth,trainingtime,fracneuron,swept_param,obsv_dstrb,'EFH');
    
    % load LDS params
    LDSparamsEM = loadEFHspikecounts(monkey,iSession,...
        binwidth,trainingtime,fracneuron,swept_param,obsv_dstrb,'LDS');
    
    % load testing spikes
    [Xtrain,Qtrain] = params.getLatents([],dataclass,...
        'sequencelength','singlesequence');
    [Rtrain,Qtrain] = params.getData(Xtrain,Qtrain);
    [Rtest,Xtest,Qtest] = params.getTestData(dataclass);
    SStot = sum((Xtest - mean(Xtest)).^2);
    
    
    keyboard
    
    % KFobs
    [rcrnt2spikeContribRatio,rcrnt2spikeUpdateContribRatio] =...
        separateSpikeAndRecurrentContributions_KFobs(Rtest,Xtest,...
        Rtrain,Xtrain,LDSparamsEM,params.Nmsperbin,SStot,TOPLOT);
    
    % KFlatent
    [rcrnt2spikeContribRatio_KF,rcrnt2spikeUpdateContribRatio_KF] =...
        separateSpikeAndRecurrentContributions_KFlatent(Rtest,Xtest,...
        Rtrain,Xtrain,LDSparamsEM,params.Nmsperbin,SStot,TOPLOT);
    
    % the rEFH
    [rcrnt2spikeContribRatio_rEFH,rcrnt2spikeUpdateContribRatio_rEFH] =...
        separateSpikeAndRecurrentContributions_rEFH(Rtest,Xtest,Qtest,...
        Rtrain,Xtrain,Qtrain,wts,params,SStot,TOPLOT);
       
    
%     %%% not ideal...
%     load('C:\#DATA\RBMish/BMI/Rsqs_BinwidthSweep_Poisson_Indy.mat','kinemat')
%     
%     fprintf('What fraction of the change in the reconstructions is due to\n')
%     fprintf(' changes in the historical input, and what fraction to changes\n')
%     fprintf(' in the spike input?  Alternatively, what is the ratio of\n')
%     fprintf(' these two contributions (historical to spike)?\n');
%     printmatJGM(...
%         [rcrnt2xhats_KF./spike2xhats_KF; rcrnt2xhats_rEFH./spike2xhats_rEFH],...
%         'foo',...
%         sprintf('%s ','KFlatent','rEFH'),...
%         sprintf('%s ',kinemat.name));
end
    
    

% plot
setColors
tikzBarGraph(...
    1:Nstates,...
    [stdZ',stdU',stdY'],...
    NaN*ones(Nstates,2,3),... % for now--put in standard error of the std later if you like
    0,...
    {'posx', 'posy', 'velx', 'vely', 'accx', 'accy'},...
    'kinematic variable',...
    sprintf('%s','units'),'',...
    repmat({'EFHclr';'EFHclr';'EFHclr'},[1,Nstates])',...
    2,0.32,'grouped',{},...
    sprintf('SLMD/%s_thing_%03d_bar',...
    'relative_contributions',binwidth));



% now plot
keyboard
plot_stuff(uz,yz,0)
plot_stuff(XhatTU,Kinnov,100)






end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function plot_stuff(running_summary,current_obsv,fig_num_base)

figure(1 + fig_num_base); clf; hold on
mu = [mean(running_summary,1); mean(current_obsv,1)]';
scatter(mu(:,1),mu(:,2));

mumu = mean(mu);
mu = mu - mumu;
[~,~,V] = svd(mu);
rr = 4; %%% fix this hack
plot([-rr*V(1,1) + mumu(1), rr*V(1,1) + mumu(2)],...
    [-rr*V(2,1) + mumu(1), rr*V(2,1) + mumu(2)])
xlabel('running summary'); ylabel('current obsv')
axis tight
hold off;


figure(2 + fig_num_base); clf; hold on;
mu2 = [mean(running_summary,2), mean(current_obsv,2)];
tvec = 500:1000; size(mu2,1);
plot(tvec,mu2(tvec,1))
plot(tvec,mu2(tvec,2))
hold off;


%%%%%
% (1) Look at FFT of fig. 2 components??

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [rcrnt2spikeContribRatio,rcrnt2spikeUpdateContribRatio] =...
    separateSpikeAndRecurrentContributions_KFobs(Rtest,Xtest,...
    Rtrain,Xtrain,LDSparamsEM,Nmsperbin,SStot,TOPLOT)

% Ns
Nsamples = size(Xtest,1);

% get params
LDSparamsOBS = fullyObservedLDSDecoder(Rtrain,Xtrain,size(Rtrain,1),...
    Rtest,Xtest,SStot);

% run the KF
LDSparamsOBS.T = Nsamples;
KFdstrbs = KalmanFilter(LDSparamsOBS,Rtest');

% get the *last* Kalman gain, after convergence
CvrnCtr = gather(KFdstrbs.INFOTU(:,:,end))\LDSparamsOBS.C';
Kfinal = CvrnCtr/(LDSparamsOBS.SigmaYX + LDSparamsOBS.C*CvrnCtr);
IminusKCfinal = eye(size(Kfinal,1)) - Kfinal*LDSparamsOBS.C;

% via the multisensory-integration interpretation of the KF
spike_contrib_to_xhat = (Rtest - LDSparamsOBS.muYX')*Kfinal';
rcrnt_contrib_to_xhat = KFdstrbs.XHATTU'*IminusKCfinal';
clear KFdstrbs
%norm(Kfinal), norm(IminusKCfinal) ?

% now summarize the contributions by their standard deviations...
rcrnt2spikeContribRatio = std(rcrnt_contrib_to_xhat)./...
    std(spike_contrib_to_xhat);
rcrnt2spikeUpdateContribRatio = std(diff(rcrnt_contrib_to_xhat))./...
    std(diff(spike_contrib_to_xhat));

end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [rcrnt2spikeContribRatio,rcrnt2spikeUpdateContribRatio] =...
    separateSpikeAndRecurrentContributions_KFlatent(Rtest,Xtest,...
    Rtrain,Xtrain,LDSparamsEM,Nmsperbin,SStot,TOPLOT)

% Ns
Nsamples = size(Xtest,1);
Nltnt = size(LDSparamsEM.C,2);
Nmspers = 1000; % fact

% run the KF
LDSparamsEM.T = Nsamples;
KFdstrbs = KalmanFilter(LDSparamsEM,Rtest');

% get the *last* Kalman gain, after convergence
CvrnCtr = gather(KFdstrbs.INFOTU(:,:,end))\LDSparamsEM.C';
Kfinal = CvrnCtr/(LDSparamsEM.SigmaYX + LDSparamsEM.C*CvrnCtr);
IminusKCfinal = eye(Nltnt) - Kfinal*LDSparamsEM.C;

% via the multisensory-integration interpretation of the KF
spike_contrib_to_latents = (Rtest - LDSparamsEM.muYX')*Kfinal';
rcrnt_contrib_to_latents = KFdstrbs.XHATTU'*IminusKCfinal';
clear KFdstrbs
%norm(Kfinal), norm(IminusKCfinal) ?

% unforunately, you don't save Bzx, so you've got to re-fit it
[~,Bzx] = latentStateLDSDecoder([],Rtrain,Xtrain,Rtest,Xtest,...
    [],[],LDSparamsEM,SStot);

% now summarize the contributions by their standard deviations...
rcrnt2spikeContribRatio = std(rcrnt_contrib_to_latents*Bzx(1:end-1,:))./...
    std(spike_contrib_to_latents*Bzx(1:end-1,:));
rcrnt2spikeUpdateContribRatio = std(diff(rcrnt_contrib_to_latents)*Bzx(1:end-1,:))./...
    std(diff(spike_contrib_to_latents)*Bzx(1:end-1,:));


% % ...or average std's over latents (which we don't care about severally)
% avg_std_spike2ltnts = mean(std(spike_contrib_to_latents));
% avg_std_rcrnt2ltnts = mean(std(rcrnt_contrib_to_latents));
% 
% 
% if TOPLOT
%     t = (0:(Nsamples-1))*Nmsperbin/Nmspers;
%     figure(1); clf; hold on;
%     iState = 3;
%     plot(t,Xtest(:,iState),'k')
%     plot(t,spike_contrib_to_latents*Bzx(1:end-1,iState))
%     plot(t,rcrnt_contrib_to_latents*Bzx(1:end-1,iState))
%     hold off;
%     legend('true','rcrnt','spike');
% end
% 

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [rcrnt2spikeContribRatio,rcrnt2spikeUpdateContribRatio] =...
    separateSpikeAndRecurrentContributions_rEFH(Rtest,Xtest,Qtest,...
    Rtrain,Xtrain,Qtrain,wts,params,SStot,TOPLOT)

% useful functions
logistic = @(XX)(1./(1 + exp(-XX)));

% Ns
Nsensory = params.numsUnits{1};
Nhids = params.numsUnits{2};
Nsamples = size(Xtest,1);

% unfortunately, you didn't save Bvx for each session
[~,~,~,~,~,Bx,~] = testEFHBMI(Rtest,Xtest,...
    Qtest,wts,params,Rtrain,Xtrain,Qtrain,SStot);

% decoding matrix
Byx = Bx(1:Nsensory,:);
Bux = Bx(Nsensory + (1:Nhids),:);
Bzx = Bx((Nsensory + Nhids) + (1:Nhids),:);
bx  = Bx(end,:);

% split up the wts into intuitive pieces
Wuz = wts{1}(1:Nhids,:);
Wyz = wts{1}((Nhids+1):(end-1),:);
bz = wts{1}(end,:);
bu = wts{2}(end,1:Nhids);
by = wts{2}(end,(Nhids+1):end);

% init/malloc
Ztest = zeros(Nsamples,Nhids,'like',Rtest);
Utest = zeros(1,Nhids,'like',Rtest);
uz = zeros(Nsamples,Nhids,'like',Rtest);
yz = Rtest*Wyz;

for iSample = 1:Nsamples
    % the recurrent "input" to the hidden layer
    uz(iSample,:) = Utest*Wuz;
    
    % compute and store the hidden activities
    Utest = logistic(uz(iSample,:) + yz(iSample,:) + bz);
    Ztest(iSample,:) = Utest;
end


% changes in output over time are small, so 
%   delta z(t) \approx f'(uz(t) + yz(t) + b)(delta uz(t) + delta yz(t))
%                   =: alpha_t*(delta uz(t) + delta yz(t))
alpha = Ztest.*(1 - Ztest);
spike_contrib_to_latents = alpha(2:end,:).*diff(yz);
rcrnt_contrib_to_latents = alpha(2:end,:).*diff(uz);

% are the changes in the output small enough to justify analysis based on
% deriatives?  Look and see:
if TOPLOT 
    figure(10);
    scatter(vect(Ztest(2:end,:)),...
        vect(spike_contrib_to_latents + rcrnt_contrib_to_latents));
    % scatter(vect(diff(Ztest)),vect(rcrnt2ltnts + spike2ltnts));
end


r = corr(vect(Ztest(2:end,:)),...
    vect(cumsum(spike_contrib_to_latents) + cumsum(rcrnt_contrib_to_latents)));
fprintf('correlation of input and output estimates is %0.3g\n', r);


r = corr(vect(diff(Ztest)),...
    vect(spike_contrib_to_latents + rcrnt_contrib_to_latents));
fprintf('correlation of input and output updates is %0.3g\n', r);


% now summarize the contributions by their standard deviations...
rcrnt2spikeContribRatio = std(cumsum(rcrnt_contrib_to_latents)*Bzx)./...
    std(cumsum(spike_contrib_to_latents)*Bzx);
rcrnt2spikeUpdateContribRatio = std(rcrnt_contrib_to_latents*Bzx)./...
    std(spike_contrib_to_latents*Bzx);


% summarize contributions of to kinematic reconstruction, xhat
%%%
% what?
%%%

% % summarize contributions of the update visible units
% XhatU = logistic(Ztest*Wuz' + bu)*Bux;
% XhatY = exp(Ztest*Wyz' + by)*Byx;
% XhatZ = Ztest*Bzx;
% stdZ = std(XhatZ);
% stdY = std(XhatY);
% stdU = std(XhatU);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function plotAutocorrelations(R,X,Nslag,Nmsperbin,fignum)

% Ns
Nmspers = 1000;
Nbins = size(X,1);
NbinsHistory = floor(Nslag*Nmspers/Nmsperbin);
samples = ((Nbins-NbinsHistory):Nbins);
t = (samples - Nbins)*Nmsperbin/Nmspers;

% plot
figure(fignum); clf; 
meanAutocorr =...
    matfun(@(iState)(mean(...
    matfun(@(iNeuron)( xcov(R(:,iNeuron), X(:,iState)) ),1:size(R,2) ),...
    2)), 1:size(X,2));
plot(t,meanAutocorr(samples,:))
legend('posx','posy','velx','vely','accx','accy')

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function XhatStaticTrain = rEFHdecodeTrainingData(wts,params,Rtrain,T,Bx)

Nunits = params.numsUnits{end};
Nlayers = length(params.numsUnits);
Nsamples = size(Rtrain,1);

[R1,Z0] = updownRDBN(Rtrain,wts,params,T);
U1 = invParamMap(Z0,wts{Nlayers}(1:end-1,1:sum(Nunits)),...
    wts{Nlayers}(end,1:sum(Nunits)),params.typeUnits{end},...
    Nunits,params);
XhatStaticTrain = [R1,U1,Z0,ones(Nsamples,1,'like',Rtrain)]*Bx;

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function Rsqs = getLagCorrelations(X,R,Nslag,Nmsperbin,titlestr)

% Ns
Nmspers = 1000;
Nbinslag = floor(Nslag*Nmspers/Nmsperbin);
Nmsperbin/Nmspers;
[Nbins,Nstates] = size(X);

% malloc
Rsqs = zeros(Nbinslag,Nstates);

% cycle through lags
for tau = 0:(Nbinslag-1)
    Mbins = Nbins - tau;
    [~,Rsqs(tau+1,:)] = linregress([ R(1:Mbins,:) ,ones(Mbins,1)], X(tau+(1:Mbins),:));
end

figure()
plot( (0:-1:(1-Nbinslag))*Nmsperbin/Nmspers, Rsqs )
ylabel('Rsq')
xlabel('time (s)')
ax = axis;
axis([ax(1), ax(2), 0, 1.0]);
title(titlestr)


end
%-------------------------------------------------------------------------%



















