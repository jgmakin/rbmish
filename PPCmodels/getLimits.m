function roboparams = getLimits(thmin,thmax,armlengths,Ndims,mods,neutralspace)
% getLimits     Get the limits of the reachable workspace
% 
% USAGES:
%   roboparams = getLimits(thmin,thmax,armlengths,params.Ndims,params.mods,params.NS)
%
% Gets the limits of the reachable workspace and assembles them as
% roboparams, along with the other facts about the "robot" arm:
%       .thmin          minimum joint angles  (shoulder,elbow)
%       .thmax          maximum joint angles  (shoulder,elbow)
%       .posmin         minimum hand position (x,y)
%       .posmax         maximum hand position (x,y)
%       .eyemin         minimum gaze angle (x,y)
%       .eyemax         maximum gaze angle (x,y)
%       .armlengths     link lengths
%       .gst            base frame
%       .w              rotation axes
%       .q              points on those axes


%-------------------------------------------------------------------------%
% Revised: 12/31/16
%   -part of massive revision to liberate generateData from case statements
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%


% store some roboparams
roboparams.thmin = thmin;
roboparams.thmax = thmax;
roboparams.armlengths = armlengths;

switch Ndims
    case 1
        roboparams.posmin = FK2link(thmax,roboparams,1)'; % -1*L1;
        roboparams.posmax = FK2link(thmin,roboparams,1)'; % 1*L1;
        if length(mods) == 3
            roboparams.eyemin = roboparams.posmin/2;
            roboparams.eyemax = roboparams.posmax/2;
            roboparams.posmin = roboparams.posmin/2;
            roboparams.posmax = roboparams.posmax/2;
        end
    case 2
        L1 = armlengths(1);                             % (upper) arm
        L2 = armlengths(2);                             % forearm
        
        % use product-of-exponentials formula
        w(:,1) = [0 0 1]';
        w(:,2) = [0 0 1]';
        q(:,1) = [0 0 0]';
        q(:,2) = [0 L1 0]';
        roboparams.gst0 = [eye(3) [0 L1+L2 0]'; 0 0 0 1];
        roboparams.w = w;
        roboparams.q = q;
        
        switch neutralspace
            case 'Joint-Angle'        
                % get position min and max
                tmp1 = FK2link([thmax(1) thmin(2)],roboparams,1)';
                tmp2 = FK2link(thmax',roboparams,1)';
                roboparams.posmin = [tmp1(1); tmp2(2)];
                tmp1 = FK2link(thmin',roboparams,1)';
                % max when second joint is as straight as possible
                thbest = thmin(2);
                % interior angle of triangle
                thint = thbest + pi/2;
                tmp2 = sqrt(L1^2 + L2^2 - 2*L1*L2*cos(thint));
                roboparams.posmax = [tmp1(1); tmp2];
                
                if length(mods) == 3
                    %%% should use roboparams rather than hard-coding 10...
                    cntr = [-10;0];                 % center of the head!
                    roboparams.eyemin = -(cntr + 10);   % your eye = -gaze
                    roboparams.eyemax = -(cntr - 10);   % your eye = -gaze
                    roboparams.posmin = roboparams.posmin + roboparams.eyemin;
                    roboparams.posmax = roboparams.posmax + roboparams.eyemax;
                    
                end
                
            case 'Hand-Position'
                
                % --- see (1) --- %
                xLowerRight = FK2link([thmin(1),thmax(2)],roboparams,1);
                gamma = cos(pi - thmin(2));
                dfr = diff(xLowerRight);
                d = 1/2*(-dfr + sqrt(2*(L1^2 + L2^2 - gamma*2*L1*L2) -...
                    sum(xLowerRight)^2));
                
                roboparams.posmin = (xLowerRight - [d 0])';
                roboparams.posmax = (xLowerRight + [0 d])';
                
                roboparams.thmin = min([IK2link(xLowerRight,roboparams,1);...
                    IK2link(xLowerRight + [-d d],roboparams,0)])';
                roboparams.thmax = max([IK2link(xLowerRight,roboparams,1);...
                    IK2link(xLowerRight + [-d d],roboparams,0)])';
                
        end
    otherwise
        error('strange number of encoded vars/neurons! -- jgm');

end

%%%% roboparams:
% posmin, posmax, thmin, thmax, eyemin, eyemax, gst0, w, q

end

% --- (1) --- %
% You want to maximize the area of a rectangle with its lower right corner
% at z, without leaving the (original) reachable workspace--i.e., the
% batman shape.  The upper left boundary (on which the point cattycorner to
% z will live) turns out to correspond to th2 at minimum (i.e., fully
% extended elbow).  We also know that the rectangle with maximum area is a
% square.  What is the length of a side of this rectangle, d?
%
% Mathematically,
%
%   IK2link_2(z + [-d,d]) = thmin(2)
%   => pi - acos(gamma) = thmin(2)
%   => gamma = cos(pi - thmin(2))
%
%   Now express d in terms of gamma:
%
%   => rsq = L1^2 + L2^2 - gamma*2*L1*L2
%   => (z1-d)^2 + (z2+d)^2 = L1^2 + L2^2 - gamma*2*L1*L2
%   => 2*d^2 + 2*(z2-z1)*d + z1^2 + z2^2 - L1^2 - L2^2 + gamma*2*L1*L2 = 0
%   => d    = [-2*(z2-z1) +/- sqrt(4*(z2-z1)^2 +...
%               8*(L1^2 + L2^2 - gamma*2*L1*L2 - z1^2 - z2^2)]/4
%           = [-(z2-z1) +/- sqrt((z2-z1)^2 +...
%               2*(L1^2 + L2^2 - gamma*2*L1*L2 - z1^2 - z2^2)]/2
%           = [-(z2-z1) +/- sqrt(-(z1+z2)^2 + 2*(L1^2 + L2^2 - gamma*2*L1*L2]/2
%
% which, somehow, magically, is 12

