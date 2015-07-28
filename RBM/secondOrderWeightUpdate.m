function [wNEW,wdotNEW] = secondOrderWeightUpdate(w,wdot,u,m,b,k,FXN,params)
% secondOrderWeightUpdate   Update neural-network weights
%
% USAGE:
%   [w,wdot] = secondOrderWeightUpdate(w,wdot,m,b,k,Ts)
%
% Rather than making weight changes ("velocities") directly proportional to
% the gradient (e.g., ML, CD-N, etc.), make the change in weight *velocity*
% (i.e., acceleration) exhibit this dependence.  In other words, make the
% gradient an "input" to a second-order discrete-time LTI system, conceived
% as a mass-spring-damper (MSD).
%
% For a longer explanation, see your labnotes.
%
%
% RELATIONSHIP TO HINTON'S PARAMETERS
%   
% If A = [1 Ts; -k/m*Ts (1-b/m*Ts)],
%
%   Ts = 1
%   m = 1/"epsilonw"
%   k = "weightcost"
%   1-b/m = "momentum"
%
% NB that changes to mass in GEH's version (e.g. the lowering of learning
% rates over time/increase in mass) do not affect the momentum, which is
% coded directly!  This implies (inter alia) that the increase in momentum
% at epoch 5 from 0.5 to 0.9 (i.e., b/m goes from 250 to 50) can be instead
% be achieved gradually by learning-rate annealing (equivalently, increase
% in mass).  [You (JGM) have traditionally multiplied the mass by sqrt(10)
% (~= 3) every ten epochs, which is fairly similar to the factor of 5
% increase implemented by the momentum change.]

%-------------------------------------------------------------------------%
% Revised: 09/30/14
%   -changed the "momentum" from (1-b/m*Ts) to (0.98 - b/m*Ts)
% Created: 09/17/14
%   by JGM
%-------------------------------------------------------------------------%

% params
Ts = params.Ts;

% weight updates
wNEW = w + Ts*wdot;
switch FXN
    case 'BP'
        % scale params by 1/alp for Bernoulli-Bernoulli units
        alp = 10; % params.g? That "translates" b/n Bernoulli and Poisson...
        t = params.t;
        wdotNEW = zeros(size(wdot),'like',wdot);
        wdotNEW(1:t,:) = -k/m*Ts*w(1:t,:) +...
            (0.98-b/m*Ts)*wdot(1:t,:) + Ts/m*u(1:t,:)*alp;
        wdotNEW(t+1:end,:) = -k/m*Ts*w(t+1:end,:) +...
            (0.98-b/m*Ts)*wdot(t+1:end,:) + Ts/m*u(t+1:end,:);
    otherwise
        wdotNEW = -k/m*Ts*w + (0.98-b/m*Ts)*wdot + Ts/m*u;
end


% *** (1) ***
% Learning in neural networks is often improved by decreasing the learning
% rate, and by increasing the "momentum" of late updates.  A single change
% can enforce both: increase the inertia of the updates.  This is performed
% outside this function, by giving the mass itself first-order linear 
% dynamics.
%
% NB that, if mass is therefore considered a third state, the total system
% is nonlinear!





















