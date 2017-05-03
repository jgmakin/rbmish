function pos = FK2link(th,roboparams,VERBOSE)
% FK2link   Forward kinematics for a 2-link robot
%
% USAGE:
% 	pos = FK2link(th,roboparams,VERBOSE)
%
% FK2link converts the joint angles in th into end-effector positions, pos.
% The oboparams structure must specify the minimum (thmax) and maximum
% (thmin) joint angles; the two axes of rotation (3x2 matrix w) and points
% through which they pass (3x2 matrix q); and the base frame (4x4 matrix
% gst0 in SE(3)).  In order to accommodate 1D kinematics, roboparams must
% also specify the length of the first link (L1).
%
% The input th can contain multiple examples of joint pairs, but must have
% size (Nexamples x Ndims)!  The output pos will have the same size.

%-------------------------------------------------------------------------%
% Revised: 12/31/16
%   -changed to run on roboparams rather than params
% Revised: 11/18/16
%   -rewrote help
% Revised: 07/02/14
%   -matricized inputs, outputs
% Revised: 11/29/10
%   -forced output to be *column* vector
% Created: 06/10/10
%   by JGM
%-------------------------------------------------------------------------%

% params
thmin = roboparams.thmin;
thmax = roboparams.thmax;
Ndims = length(thmin);

% check for out-of-range angles
if VERBOSE
    outofrange(th,thmin,thmax,'Joint-Angle',1);
end
    
% 1- or 2-D kinematics
if Ndims == 1
    pos = cos(th)*roboparams.armlengths(1);
else
    % omega = axis of rotation, q = pt. on that axis
    w1 = roboparams.w(:,1); q1 = roboparams.q(:,1);
    w2 = roboparams.w(:,2); q2 = roboparams.q(:,2);
    gst00 = roboparams.gst0;  
    % maps pts in S frame into the T frame @ th = 0
    
    % compute screw motions
    screw1 = screw(w1,cross(q1,w1),th(:,1));
    screw2 = screw(w2,cross(q2,w2),th(:,2));
    pos02 = tensorOp(screw2,gst00(:,end,ones(size(th,1),1)));
    pos12 = tensorOp(screw1,pos02);
    
    % you only need the first two components
    pos = permute(pos12(1:2,:,:),[3,1,2]);
    
end

end

%%% gst = screw1*screw2*gst0;
%%% pos(1,1) = gst(1,end);
%%% pos(2,1) = gst(2,end);