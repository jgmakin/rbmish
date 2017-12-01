function [beta_best,RsqWF] = WFDecoder(Rtrain,Xtrain,T,N_ms_per_bin,...
    Rtest,Xtest)
% WFDecoder
%
% USAGE....
%
% NB that Nhistorybins = 1 means just using the current sample

%-------------------------------------------------------------------------%
% Cribbed: 09/29/17
%   -from filtersForNeuralData.m (JGM)
% Created: ~04/xx/16
%   by JGM
%-------------------------------------------------------------------------%


% Ns
N_samples = size(Rtrain,1);
N_seq = floor(N_samples/T);
N_ms_per_s = 1000;

% reshape data into trajectories
Xtrain = shortdata(N_seq,3,Xtrain);
Rtrain = shortdata(N_seq,3,Rtrain);

% useful params to be set
frac_train = 0.8; %%% make into an argument?
N_bins_step = 1;
N_sec_history_max = 0.4;
%%%lambda_min = 800;
lambda_min = 0;
lambda_max = 4000;
lambda_step = 800;
TOPLOT = 0;

% get ranges of lambdas and history lengths (in number of bins)
N_history_bins_max = round(N_sec_history_max/(N_ms_per_bin/N_ms_per_s));
all_lambdas = lambda_min:lambda_step:lambda_max;
all_N_history_bins = 1:N_bins_step:N_history_bins_max;
RsqWF_heldout = zeros(size(Xtrain,2), length(all_N_history_bins),...
    length(all_lambdas), 'like', Xtrain);

% to get a baseline for performance, use the mean
%%% You might think you should acquire the means on the first frac_train of
%%% the data, then evaluate the corresponding best_vel_Rsq on the remaining
%%% training data--but that's wrong: you *know* what using the means should
%%% give you--R^2=0 for all vars--so you don't have to evaluate anything.
%%% Moreover, an evaluation on just 20% of the training data will be noisy
%%% anyway; *and* holding no data out gives a better estimate of the means.
beta_best = zeros(size(Rtest,2), size(Xtest,2), 'like', Rtest);
beta_best = [beta_best; mean(mean(Xtrain,3),1)];
best_vel_Rsq = 0;
lambda_best = 0;


fprintf('Sweeping hyperparams for Wiener filter...\n');
for i_lambda = 1:length(all_lambdas)
    for i_NhistoryBins = 1:length(all_N_history_bins)
        [this_beta, this_Rsq] = LearnWienerFilter(Rtrain,Xtrain,...
            all_N_history_bins(i_NhistoryBins), all_lambdas(i_lambda),...
            frac_train);
        
        % maximize *velocity* decoding
        this_vel_Rsq = mean(this_Rsq(3:4));
        if this_vel_Rsq > best_vel_Rsq
            best_vel_Rsq = this_vel_Rsq;
            beta_best = this_beta;
            lambda_best = all_lambdas(i_lambda);
        end
        
        % store *all* Rsqs
        RsqWF_heldout(:,i_NhistoryBins,i_lambda) = this_Rsq';
        fprintf('.')
    end
end
fprintf('\n')

if TOPLOT
    fprintf('best regularization is %f\n',lambda_best);
    vars = {'x', 'y', 'xdot', 'ydot', 'xddot', 'yddot'};
    for iVar = 1:size(RsqWF_heldout,1)
        figure(iVar)
        plot(all_N_history_bins*N_ms_per_bin,squeeze(RsqWF_heldout(iVar,:,:)));
        title(vars{iVar});
        ax = axis; axis([ax(1:3), 1.0])
        ylabel('Rsq')
        xlabel('amount of history (ms)')
        legend(arrayfun(@(i)(['\lambda = ',num2str(all_lambdas(i))]),...
            1:length(all_lambdas), 'UniformOutput', false),...
            'location', 'northwest')
    end
end

% % for each kinematic variable, get best lambda, Nhistorybins
% Nvars = size(RsqWF,1);
% history_argmax = zeros(Nvars,1);
% lambda_argmax = zeros(Nvars,1);
% for iVar = 1:Nvars
%     these_Rsqs = RsqWF(iVar,:,:);
%     [~,Rsq_max_ind] = max(these_Rsqs(:)); 
%     [Rsq_history_max_ind,Rsq_lambda_max_ind] =...
%         ind2sub(size(squeeze(these_Rsqs)),Rsq_max_ind);
%     history_argmax(iVar) = all_Nhistorybins(Rsq_history_max_ind);
%     lambda_argmax(iVar) = all_lambdas(Rsq_lambda_max_ind);
% end

clear Xtrain Rtrain

% now compute Rsqs on the *test* data using the best betas
Nhistorybins_best = (size(beta_best,1)-1)/size(Rtest,2);
[Rtest, Xtest] = shortdata(N_seq,3,Rtest,Xtest);
[Rtest, Xtest] = historicizeSampleData(Rtest,Xtest,Nhistorybins_best);
Xres = Xtest - Rtest*beta_best;
SStot = sum((Xtest - mean(Xtest,1)).^2);
SSerr = sum(Xres.^2,1);
RsqWF = 1 - SSerr./SStot;



%%% fix this to accommodate Xhats w/fewer than the full set of samples
%if ~isempty(tvec), PosVelPlot(tvec,aa*beta_best,0); end
%%% 


end
