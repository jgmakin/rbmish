function [wNEW,wdotNEW] = secondOrderWeightUpdate(w,wdot,u,m,b,k,nums,params)
% secondOrderWeightUpdate   Update neural-network weights
%
% USAGE:
%   [w,wdot] = secondOrderWeightUpdate(w,wdot,u,m,b,k,numsUnits,params);
%
% Rather than making weight changes ("velocities") directly proportional to
% the gradient (e.g., ML, CD-N, etc.), make the change in weight *velocity*
% (i.e., acceleration) exhibit this dependence.  In other words, make the
% gradient, u, an "input" to a second-order discrete-time LTI system, 
% conceived as a mass-spring-damper (MSD).
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
%   1-b/m*Ts = "momentum"
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
%%% wNEW = w + Ts*wdot;
if length(unique(m)) > 1
    endinds = cumsum(nums);
    startinds = [1, endinds(1:end-1)+1];
    wdotNEW = zeros(size(wdot),'like',wdot);
    for iGrp = 1:length(m)
        inds = startinds(iGrp):endinds(iGrp);
        mm = m(iGrp)/Ts;                % normalize units so that Ts = 1
        wdotNEW(inds,:) = -k(iGrp)/mm*w(inds,:) +...
            (0.98-b(iGrp)/mm)*wdot(inds,:) + u(inds,:)/mm;
    end
else
    mm = m(1)/Ts;                       % normalize units so that Ts = 1
    wdotNEW = (-k(1)/mm)*w + (0.98-b(1)/mm)*wdot + u/mm;
end
wNEW = w + Ts*wdotNEW;

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

% *** (2) ***
% In the recurrent models, changing the weights changes changes the "data."
% Thus, gradients measured under one set of weights should have their
% effect *on that set of weights*.  That requires updating the position
% immediately, rather than waiting one iteration.  So now technically it
% doesn't look precisely like a 2nd-order dynamical system, at least not
% with the positions and velocities grouped this way.
















