clear; clc; close all

% load
% load results/wts15smpllrn.mat
% load results/NN15smplLrn.mat
% load results/new/wtsBigSpace
load results/finalwts/Std050.mat

RESTART = 1;
NSIND = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELETE ME? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D0, S0] = generateData(1000,params);
SINSMerr = covInCalc(D0,S0,params);
ErrCovIn = SINSMerr{1}.cov + SINSMerr{2}.cov;
iDirection = 0; nDirections = 8;
phi = 2*pi*iDirection/nDirections;
R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
shft = sqrtm(ErrCovIn)*R*[2.5; 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shft = [0;0];


%%
% calculate
if RESTART
    i = 1;
    ErrorStatsArray = [];
    for pgain = 6:3:24 % 5:2:25
        vgain = 30 - pgain;
        
        % test
        params.swing = 0;
        ErrorStats = test(wts,params,'propbias',shft,...
            'propgain',pgain,'visgain',vgain,...
            'propagation','Nsamples','numsamples',15);
        
        gainratios(i) = pgain/vgain;
        ErrorStatsArray = [ErrorStatsArray; ErrorStats];
        
        i=i+1;
    end  
    save results/gainmodstats0 ErrorStatsArray gainratios params
else
    load results/new/gainmodstats
end


%%
% plot
close all


% figure(); hold on;
LL = length(gainratios);
hsp = subplot(2,(LL+1)/2,(LL+1)/2+1);
eRBM = ErrorStatsArray{(LL+1)/2,3};
eOPT = ErrorStatsArray{(LL+1)/2,4};
hold on;
h = error_ellipse(eRBM{NSIND}.cov/40000,eRBM{NSIND}.mu,'conf',.95);
set(h,'LineWidth',1,'Color','g');
h = error_ellipse(eOPT{NSIND}.cov/40000,eOPT{NSIND}.mu,'conf',.95);
set(h,'LineWidth',1,'Color','k');
axis equal
p = get(hsp,'pos');
p(3) = p(3) + 0.03;
set(hsp,'pos',p);


vec = [(LL+1)/2:-1:1,(LL+1)/2+2:1:LL+1];
for i = 1:LL
    eRBM = ErrorStatsArray{i,3};
    eOPT = ErrorStatsArray{i,4};
  
    hsp = subplot(2,6,vec(i));
    hold on;
    h = error_ellipse(eRBM{NSIND}.cov/40000,eRBM{NSIND}.mu,'conf',.95);
    set(h,'LineWidth',1,'Color','g');% [i/LL 1-i/LL i/LL]);
    h = error_ellipse(eOPT{NSIND}.cov/40000,eOPT{NSIND}.mu,'conf',.95);
    set(h,'LineWidth',1,'Color','k'); % [i/LL 1-i/LL i/LL]);
    axis equal
    p = get(hsp,'pos');
    p(3) = p(3) + 0.03;
    set(hsp,'pos',p);
    hold off;
    
    biasRBM(:,i) = eRBM{NSIND}.mu - eOPT{NSIND}.mu;
    ecovRBM(:,:,i) = eRBM{NSIND}.cov;
    detRBM(i) = det(eRBM{NSIND}.cov);
    detOPT(i) = det(eOPT{NSIND}.cov);
    
    [~, D] = eig(eOPT{NSIND}.cov\eRBM{NSIND}.cov);
    covdist(i) = norm(log(diag(D)));
    
end
% legend('5:25','7:23','9:21','11:19','13:17','15:15','17:13','19:11','21:9','23:7','25:5');
% for i = 1:LL
%     eRBM = ErrorStatsArray{i,3};
%     scatter(eRBM{NSIND}.mu(1),eRBM{NSIND}.mu(2),'kx');
%     MUS(:,i) = eRBM{NSIND}.mu;
% end
% scatter(0,0,'rx')
% title('Biases')
% plot(MUS(1,:),MUS(2,:),'k')
% hold off;
% figure; hold on;
% for i = 1:LL
%     
%     eOPT = ErrorStatsArray{i,4};
%    
%     
% end

%%
RBMbmag = sqrt(diag(biasRBM'*biasRBM));
figure; hold on;
plot(log(gainratios),RBMbmag)
scatter(log(gainratios),RBMbmag)
xlabel('log(gp/gv)')
ylabel('|bias|');
title('Network Bias')

figure; hold on;
plot(log(gainratios),(detRBM-detOPT)./detOPT);
scatter(log(gainratios),(detRBM-detOPT)./detOPT);
xlabel('log(gp/gv)')
ylabel('(|RBM|-|OPT|)/|OPT|');
title('covs'); hold off;

figure; hold on;
plot(log(gainratios),covdist);
scatter(log(gainratios),covdist);
xlabel('log(gp/gv)')
ylabel('cov discrepancy');
title('Network Covariance Discrepancy'); hold off;


% confidence limits
% n = 40000;                                      %%% hard-coded
% std = chol(vecCov);
% c = 1.96;                                       % from the CDF for a normal
% figure; hold on; 
% h = errorbar(vecMu,c*std/sqrt(n),'xr','LineWidth',1);
% plot(zeros(1,length(vecMu)),'k');          
% scatter(1:length(vecCov),vecCov/5,'k')
% hold off;
% 



%%
% plot
% close all
clc
setColors

close all
f1 = figure(1);
pstn = get(f1,'pos');
set(f1,'pos',pstn + [-400 0 600 0]);

f2 = figure(2);
pstn = get(f2,'pos');
set(f2,'pos',pstn + [-400 0 600 0]);


for i = 1:length(gainratios)
    eINP = ErrorStatsArray{i,1};
    eRBM = ErrorStatsArray{i,3};
    eOPT = ErrorStatsArray{i,4};
  
    figure(1)
    hsp = subplot(2,6,i); % length(gainratios),i);
    hold on;
    h = error_ellipse(eINP{1}.cov,eINP{1}.mu,'conf',.95);
    set(h,'LineWidth',1,'Color',VIScolor);
    h = error_ellipse(eINP{NSIND}.cov,eINP{NSIND}.mu,'conf',.95);
    set(h,'LineWidth',1,'Color',PROPcolor);
    h = error_ellipse(eRBM{NSIND}.cov,eRBM{NSIND}.mu,'conf',.95);
    set(h,'LineWidth',1,'Color',EFHcolor);
    h = error_ellipse(eOPT{NSIND}.cov,eOPT{NSIND}.mu,'conf',.95);
    set(h,'LineWidth',1,'Color',OPTcolor);
    hold off;
    % axis([-0.11 0.11 -0.11 0.11])
    % axis equal
    
    
    axis([-0.11 0.18 -0.18 0.11])
    % axis equal

    figure(2)
    hsp = subplot(2,6,i);
    hold on;
    h = error_ellipse(eRBM{NSIND}.cov/40000,eRBM{NSIND}.mu,'conf',.95);
    set(h,'LineWidth',1,'Color','g');% [i/LL 1-i/LL i/LL]);
    h = error_ellipse(eOPT{NSIND}.cov/40000,eOPT{NSIND}.mu,'conf',.95);
    set(h,'LineWidth',1,'Color','k'); % [i/LL 1-i/LL i/LL]);
    hold off;
    % axis([-14e-4 3.1e-4 -5e-4 4e-4])
    axis equal
    % axis([0.085 0.1214 -0.1304 -0.097]);
    
    
end
