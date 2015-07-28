function pos = FK2link(th,params,VERBOSE)
% FK2link forward kinematics for a 2-link robot
%   USAGE:
%       pos = FK2link(th,params,VERBOSE)
%
%   FK2link generates random joint angles ANGLE and corresponding end
%   effector position POS for a two-link arm.  The joint limits, link
%   lengths, axis directions and locations, and initial position are set in
%   the first few lines.
%
%   USAGE:
%       pos = FK2link(th,params,VERBOSE)
%
%   NB: INPUT TH *MUST* BE IN LONGDATA FORMAT (Nexamples x Ndims)!

%-------------------------------------------------------------------------%
% Revised: 07/02/14
%   -matricized inputs, outputs
% Revised: 11/29/10
%   -forced output to be *column* vector
% Created: 06/10/10
%   by JGM
%-------------------------------------------------------------------------%

% check for out-of-range angles
if VERBOSE
    outofrange(th,params.thmin,params.thmax,'Joint-Angle',1,params);
end
    
% 1- or 2-D kinematics
if params.Ndims == 1
    pos = cos(th)*params.L1;
else
    % init
    w1 = params.w(:,1); q1 = params.q(:,1);     % omega = axis of rotation
    w2 = params.w(:,2); q2 = params.q(:,2);     % q = pt. on that axis
    gst00 = params.gst0;                         % maps pts in S frame into
                                                %  the T frame @ th = 0
    % compute screw motions
    screw1 = screw(w1,cross(q1,w1),th(:,1));
    screw2 = screw(w2,cross(q2,w2),th(:,2));
    gstpt00 = gst00(:,end);
    pos02 = tensorOp(screw2,gstpt00(:,ones(size(th,1),1)));
    pos12 = tensorOp(screw1,pos02);
    
    % you only need the first two components
    pos = pos12(1:2,:)';

end

end

%%% gst = screw1*screw2*gst0;
%%% pos(1,1) = gst(1,end);
%%% pos(2,1) = gst(2,end);