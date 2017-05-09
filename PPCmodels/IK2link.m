function th = IK2link(pos,roboparams,VERBOSE)
% IK2link Inverse kinematics for a 2-link robot
%   USAGE:
%       th = IK2link(pos,roboparams,VERBOSE)
%
%   IK2link is the complementary file to FK2link.

%-------------------------------------------------------------------------%
% Revised: 07/03/14
%   -matricized input pos, output th
% Revised: 01/07/14
%   -changed 'pos' to 'Hand-Position' and 'th' to 'Joint-Angle'
% Created: 08/18/10
%   by JGM
%-------------------------------------------------------------------------%


% Ns
Ndims = length(roboparams.thmin);

% check for out-of-range positions
if VERBOSE
    % if you sampled uniformly in prop, these won't be used
    inmin = roboparams.posmin;
    inmax = roboparams.posmax;
    if isfield(roboparams,'eyemin')
        inmin = inmin + roboparams.eyemin;
        inmax = inmax + roboparams.eyemax;
    end
    outofrange(pos,inmin,inmax,'Hand-Position',1);
end


if Ndims == 1
    th = acos(pos/roboparams.armlengths(1));  
else
    % init
    L1 = roboparams.armlengths(1);
    L2 = roboparams.armlengths(2);
    
    % intermediate params
    phi = atan2(-pos(:,1),pos(:,2));
    
    rsq = pos(:,1).^2 + pos(:,2).^2;
    gamma = (L1^2 + L2^2 - rsq)/(2*L1*L2);
    delta = (rsq + L1^2 - L2^2)./(2*L1*sqrt(rsq));
    
    % cap these, lest you get imaginary numbers somewhere
    gamma = sign(gamma).*(abs(gamma)>=1) + (abs(gamma)<1).*gamma;
    delta = sign(delta).*(abs(delta)>=1) + (abs(delta)<1).*delta;
    
    alpha = acos(gamma);
    beta = acos(delta);
    

    if ~all(isreal(beta)) || ~all(isreal(alpha))
        fprintf('warning: imaginary numbers in inverse kinematics\n');
    end

    % theta
    th(:,1) = phi - beta;
    th(:,2) = pi - alpha;
    
end
    
% % check for out-of-range angles
% if VERBOSE
%     outofrange(th,params.thmin,params.thmax,'Joint-Angle',1);
% end


end