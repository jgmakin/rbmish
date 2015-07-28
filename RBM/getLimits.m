function params = getLimits(params)


switch params.Ndims
    case 1
        params.posmin = FK2link(params.thmax,params,1)'; % -1*L1;
        params.posmax = FK2link(params.thmin,params,1)'; % 1*L1;
        if params.Nmods == 3
            params.eyemin = params.posmin/2;
            params.eyemax = params.posmax/2;
            params.posmin = params.posmin/2;
            params.posmax = params.posmax/2;
        end
    case 2
        L1 = params.L1;                                 % (upper) arm
        L2 = params.L2;                                 % forearm ;-|
        thmin = params.thmin;
        thmax = params.thmax;
        
        % use product-of-exponentials formula
        w(:,1) = [0 0 1]';
        w(:,2) = [0 0 1]';
        q(:,1) = [0 0 0]';
        q(:,2) = [0 params.L1 0]';
        params.gst0 = [eye(3) [0 L1+L2 0]'; 0 0 0 1];
        params.w = w;
        params.q = q;
        
        switch params.NS
            case 'Joint-Angle'        
                % get position min and max
                tmp1 = FK2link([thmax(1) thmin(2)],params,1)';
                tmp2 = FK2link(thmax',params,1)';
                params.posmin = [tmp1(1); tmp2(2)];
                tmp1 = FK2link(thmin',params,1)';
                % max when second joint is as straight as possible
                thbest = thmin(2);
                % interior angle of triangle
                thint = thbest + pi/2;
                tmp2 = sqrt(L1^2 + L2^2 - 2*L1*L2*cos(thint));
                params.posmax = [tmp1(1); tmp2];
                
                if params.Nmods == 3
                    %%% should use params rather than hard-coding 10...
                    cntr = [-10;0];                 % center of the head!
                    params.eyemin = -(cntr + 10);   % your eye = -gaze
                    params.eyemax = -(cntr - 10);   % your eye = -gaze
                    params.posmin = params.posmin + params.eyemin;
                    params.posmax = params.posmax + params.eyemax;
                    
                end
                
            case 'Hand-Position'
                
                % --- see (1) --- %
                xLowerRight = FK2link([thmin(1),thmax(2)],params,1);
                gamma = cos(pi - thmin(2));
                dfr = diff(xLowerRight);
                d = 1/2*(-dfr + sqrt(2*(L1^2 + L2^2 - gamma*2*L1*L2) -...
                    sum(xLowerRight)^2));
                
                params.posmin = (xLowerRight - [d 0])';
                params.posmax = (xLowerRight + [0 d])';
                
                params.thmin = min([IK2link(xLowerRight,params,1);...
                    IK2link(xLowerRight + [-d d],params,0)])';
                params.thmax = max([IK2link(xLowerRight,params,1);...
                    IK2link(xLowerRight + [-d d],params,0)])';
                
        end
    otherwise
        error('strange number of encoded vars/neurons! -- jgm');

end



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

