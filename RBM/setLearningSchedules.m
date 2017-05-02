function params = setLearningSchedules(mw0,mvb0,mhb0,GROWTH,params,mrw0,mrb0)
% setLearningSchedules  Set EFH learning-rate annealing schedules
%
% USAGES:
%   params = setLearningSchedules(500,120,120,'exp',params)
%   params = setLearningSchedules(200,200,400,'hyperbolic',params);
%   params = setLearningSchedules(500,120,120,'logistic',params);
%
% Consider the inverse of the learning rate, which we can treat as a mass.
% Then the "annealing" of the learning rate is equivalent to the mass grow-
% ing according to some dynamics.  In particular, "momentum" and "weight
% decay" of the NN parameter updates can be interpreted as facts about a
% second-order dynamical system, whose input is the gradient and whose
% output is the parameter update--see secondOrderWeightUpdate.m for
% details.  Then the mass itself grows across parameter updates.
%
% NB: if the damping coefficient, b, is set to half the initial mass (as is
% usually often down below), then oscillations can be prevented by setting
% k < m/16---at least for the continuous-time version.  In discrete time
% your margin gets slightly bigger.
%
% For Poisson units (with gains around 10?), a good choice for mw0 is 500, 
% and for the (Bernoulli) biases, 120.  For Bernoulli units, mw0 can be set
% to ~50.  Presumably this ratio is controlled by g.

%-------------------------------------------------------------------------%
% Revised: 12/26/16
%   -removed default scaling for "extra" typeUnits (assumed to be recurrent
%   units), in favor of explicitly setting with vector-valued masses.
% Created: 08/25/16 (happy anniversary, TAM&DMM)
%   by JGM
%-------------------------------------------------------------------------%


% scale down the recurrent mass (and damping)
if length(params.numsUnits{1})~=length(mw0)||...
        length(params.numsUnits{1})~=length(mvb0)||...
        length(params.numsUnits{2})~=length(mhb0)
    error(['the number of typeUnits doesn''t match the number of',...
        'initial learning ''masses'' -- jgm']);
end

% "sampling interval" in "seconds"/weight update (so m,b,k are in SI units)
Ts = 0.01;                                  

switch GROWTH
    
    case 'exp'
        
        % mass grows exponentially (learning rate decays exponentially)
        amass = 1.10;                           % State-transition "matrix"
        expGrowth = @(x0,ii)(amass^ii*x0*Ts^2);
        
        % masses grows exponentially
        params.mw  = @(iEp)(expGrowth(mw0,iEp));
        params.mvb = @(iEp)(expGrowth(mvb0,iEp));
        params.mhb = @(iEp)(expGrowth(mhb0,iEp));
        
        % damping coefficients are fixed at m0/2
        bw = mw0/2;     bvb = mvb0/2;   bhb = mhb0/2;
        params.bw  = @(iEp)(bw*Ts);
        params.bvb = @(iEp)(bvb*Ts);
        params.bhb = @(iEp)(bhb*Ts);
        
        % spring constants are "small"  -- function for generality
        params.kw  = @(iEp)(mw0/500000);
        params.kvb = @(iEp)(mvb0/500000);
        params.khb = @(iEp)(mhb0/500000);
        %%% originally 0.001 for all three--these haven't been checked
        
        if nargin > 5
            params.mrw = @(iEp)(expGrowth(mrw0,iEp));
            params.mrb = @(iEp)(expGrowth(mrb0,iEp));
            brw = mrw0/2;     brb = mrb0/2;
            params.brw = @(iEp)(brw*Ts);
            params.brb = @(iEp)(brb*Ts);
            params.krw = @(iEp)(mrw0/500000);
            params.krb = @(iEp)(mrb0/500000);
        end
        
        
    case 'logistic'
        
        % mass grows logistically
        logisticGrowth = @(x0,ii)(1000/(1 + exp(-ii/8 + 7.5)) + 0.5)*x0*Ts^2;
        
        % mass grows according to the logistic function
        params.mw  = @(iEp)(logisticGrowth(mw0,iEp));
        params.mvb = @(iEp)(logisticGrowth(mvb0,iEp));
        params.mhb = @(iEp)(logisticGrowth(mhb0,iEp));
        
        % damping coefficients are fixed at m0/2
        bw = mw0/2;     bvb = mvb0/2;   bhb = mhb0/2;
        params.bw  = @(iEp)(bw*Ts);
        params.bvb = @(iEp)(bvb*Ts);
        params.bhb = @(iEp)(bhb*Ts);
        
        % spring constants are "small"
        params.kw  = @(iEp)(0.001);
        params.kvb = @(iEp)(0.001);
        params.khb = @(iEp)(0.001);
        
        if nargin > 5
            params.mrw = @(iEp)(logisticGrowth(mrw0,iEp));
            params.mrb = @(iEp)(logisticGrowth(mrb0,iEp));
            brw = mrw0/2;     brb = mrb0/2;
            params.brw = @(iEp)(brw*Ts);
            params.brb = @(iEp)(brb*Ts);
            params.krw = @(iEp)(0.001);
            params.krb = @(iEp)(0.001);
        end
        
        
    case 'hyperbolic'
        
        % mass grows hyperbolically (learning rates decays linearly)
        hyperbolicGrowth = @(x0,ii)(x0*Ts^2/(1-(ii-1)/params.NepochsMax));
        
        % mass grows according to a hyperbola
        params.mw = @(iEp)(hyperbolicGrowth(mw0,iEp));
        params.mvb = @(iEp)(hyperbolicGrowth(mvb0,iEp));
        params.mhb = @(iEp)(hyperbolicGrowth(mhb0,iEp));
        
        % damping coefficients also follow a hyperbola, so the *momentum* 
        % is fixed---at 1-b/m*Ts = 1 - 0.08 = 0.92.
        params.bw = @(iEp)(hyperbolicGrowth(0.08*mw0/Ts,iEp));
        params.bvb = @(iEp)(hyperbolicGrowth(0.08*mvb0/Ts,iEp));
        params.bhb = @(iEp)(hyperbolicGrowth(0.08*mhb0/Ts,iEp));
        
        % no weight decay
        params.kw = @(iEp)(0);
        params.kvb = @(iEp)(0);
        params.khb = @(iEp)(0);
        
        % learning rates for recurrent units
        if nargin > 5
            params.mrw = @(iEp)(hyperbolicGrowth(mrw0,iEp));
            params.mrb = @(iEp)(hyperbolicGrowth(mrb0,iEp));
            params.brw = @(iEp)(hyperbolicGrowth(0.08*mrw0/Ts,iEp));
            params.brb = @(iEp)(hyperbolicGrowth(0.08*mrb0/Ts,iEp));
            params.krw = @(iEp)(0);
            params.krb = @(iEp)(0);
        end
        
    otherwise
        error('unexpected learning schedule! -- jgm');
end

params.Ts = Ts;
fprintf('Using the %s learning-rate adjustment scheme!!\n',GROWTH);

end
