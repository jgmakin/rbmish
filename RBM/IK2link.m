function th = IK2link(pos,params,VERBOSE)
% IK2link Inverse kinematics for a 2-link robot
%   USAGE:
%       th = IK2link(pos,params,VERBOSE)
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

% init
L1 = params.L1;

% check for out-of-range positions
if strcmp(params.MODEL,'1Daddition') || strcmp(params.MODEL,'2Daddition')
    inmin = params.posmin + params.eyemin;
    inmax = params.posmax + params.eyemax;
else
    inmin = params.posmin;      % if you sampled uniformly in prop space,
    inmax = params.posmax;      %  these won't be used
end
if VERBOSE
    outofrange(pos,inmin,inmax,'Hand-Position',1,params);
end


if params.Ndims == 1
    th = acos(pos/L1);  
else
    % init
    L2 = params.L2;
    
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
%     outofrange(th,params.thmin,params.thmax,'Joint-Angle',1,params);
% end


end