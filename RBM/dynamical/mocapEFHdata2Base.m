function Motion = mocapEFHdata2Base(skel,data)
% mocapEFHdata2Base Convert motion-capture EFH data into body-coord version
% 
% USAGE:
%   newdata = mocapEFHdata2Base(skel,Q)
%
% where Q is a structure with fields R, seqindex, offsets, data_mean, and
% data_std.  See mocapBase2EFHdata.m.

%-------------------------------------------------------------------------%
% Revised: 11/21/16
%   -functionized, tensorized, rationalized
%   by JGM
%
%-------------------------------------------------------------------------%
%
% Version 1.000
%
% Code provided by Graham Taylor, Geoff Hinton and Sam Roweis
%
% For more information, see:
%     http://www.cs.toronto.edu/~gwtaylor/publications/nips2006mhmublv
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with
% a note saying that the original programs are available from our
% web page.
% The programs and documents are distributed without any warranty, express or
% implied.  As the programs were written for research purposes only, they have
% not been tested to the degree that would be advisable in any important
% application.  All use of these programs is entirely at the user's own risk.
%
% This program postprocesses data.
% It is the first stage of preprocessing (essentially reverses
% preprocess2.m)
%
% We support two types of skeletons:
%  1) Those built from the CMU database (acclaim)
%     http://mocap.cs.cmu.edu/
%  2) Those built from data from Eugene Hsu (mit)
%     http://people.csail.mit.edu/ehsu/work/sig05stf/
%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
% [[Sutskever? notes]]
% Note our convention (because I started this way with the CMU code)
% x- ground plane axis (looking at screen, going to the right)
% z- ground plane axis (looking at screen, going "in")
% y- vertical axis
% For MIT data, subject at 0 position is facing positive x axis
% channel 1 - x component (exponential map)
% channel 2 - z component (exponential map)
% channel 3 - y (vertical) component (exponential map)
% Note this is the way that MATLAB expects it -- for the CMU data, we need
%   to reverse channels 2,3 when we plot because the ordering is (x,y,z)
%-------------------------------------------------------------------------%


for iSeq = 1:length(data.seqindex)
    
    % rescale into the original range, ...
    dataBC = EFHdata2bodyCoords(skel,data.R(data.seqindex{iSeq},:),data);
    
    % convert velocities to positions (integrate)
    [dataBase, yaw] = vel2pos(skel,dataBC);
    
    % convert from body-centred orientations to exponential maps [??]
    Motion{iSeq} = BC2expCoords(skel,dataBase,yaw);
    
end


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function dataBC = EFHdata2bodyCoords(skel,visible,data)

% scale back to original
dataBC = repmat(data.sig, size(visible,1),1) .* ...
    visible + repmat(data.mu, size(visible,1),1);

switch skel.type
    case 'acclaim'  % CMU-style data
        %No offsets, but several dimensions are constant
        
        dataBC = ...
            [ dataBC(:,1) dataBC(:,2) dataBC(:,3) ...%x,z,y rot
            dataBC(:,4) dataBC(:,5) dataBC(:,6)... %x,y,z pos
            zeros(size(dataBC,1),3) dataBC(:,7:10) ...
            zeros(size(dataBC,1),2) dataBC(:,11:14) ...
            zeros(size(dataBC,1),5) dataBC(:,15:18) ...
            zeros(size(dataBC,1),2) dataBC(:,19:22) ... %end right foot
            zeros(size(dataBC,1),2) dataBC(:,23:40) ... %end head
            zeros(size(dataBC,1),3) dataBC(:,41:44) ... %end lradius
            zeros(size(dataBC,1),3) dataBC(:,45) ... %end lwrist
            zeros(size(dataBC,1),1) dataBC(:,46:48) ...
            zeros(size(dataBC,1),3) dataBC(:,49:51) ... %end lthumb
            zeros(size(dataBC,1),3) dataBC(:,52:55) ... %end rwrist
            zeros(size(dataBC,1),3) dataBC(:,56) ... %end rwrist
            zeros(size(dataBC,1),1) dataBC(:,57:59) ...
            zeros(size(dataBC,1),3) dataBC(:,60:62)];
        
    case 'mit'
        %MIT-style data
        
        dataBC = ...
            [ dataBC(:,1) dataBC(:,2) dataBC(:,3) ...%x,z,y rot
            dataBC(:,4) dataBC(:,5) dataBC(:,6)... %x,y,z pos
            dataBC(:,7:9) repmat(data.offsets(1,:),size(dataBC,1),1) ...
            zeros(size(dataBC,1),1) dataBC(:,10) ...
            zeros(size(dataBC,1),1) repmat(data.offsets(2,:),size(dataBC,1),1) ...
            dataBC(:,11:13) repmat(data.offsets(3,:),size(dataBC,1),1) ...
            zeros(size(dataBC,1),1) ...
            dataBC(:,14) zeros(size(dataBC,1),1) ...
            repmat(data.offsets(4,:),size(dataBC,1),1) ...
            dataBC(:,15:17) repmat(data.offsets(5,:),size(dataBC,1),1) ...
            zeros(size(dataBC,1),1) ...
            dataBC(:,18) zeros(size(dataBC,1),1) ...
            repmat(data.offsets(6,:),size(dataBC,1),1) ...
            dataBC(:,19:21) repmat(data.offsets(7,:),size(dataBC,1),1) ...
            zeros(size(dataBC,1),1) ...
            dataBC(:,22) zeros(size(dataBC,1),1) ...
            repmat(data.offsets(8,:),size(dataBC,1),1) ...
            dataBC(:,23:25) repmat(data.offsets(9,:),size(dataBC,1),1) ...
            dataBC(:,26:28) repmat(data.offsets(10,:),size(dataBC,1),1) ...
            dataBC(:,29:31) repmat(data.offsets(11,:),size(dataBC,1),1) ...
            dataBC(:,32:34) repmat(data.offsets(12,:),size(dataBC,1),1) ...
            dataBC(:,35:37) repmat(data.offsets(13,:),size(dataBC,1),1) ...
            dataBC(:,38:40) repmat(data.offsets(14,:),size(dataBC,1),1) ...
            dataBC(:,41:43) repmat(data.offsets(15,:),size(dataBC,1),1) ...
            dataBC(:,44:46) repmat(data.offsets(16,:),size(dataBC,1),1) ...
            dataBC(:,47:49) repmat(data.offsets(17,:),size(dataBC,1),1)];
        
    otherwise
        error('Unknown skeleton type');
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [dataBase, yaw] = vel2pos(skel,dataBC)

% init: most dimensions of dataBC won't be changed
dataBase = dataBC;

% extract the velocities
vertrotdelta = dataBC(:,3);
grounddelta = dataBC(:,4:5);

% the original dimensions are different for CMU and MIT
switch skel.type
    case 'acclaim'  % CMU DATA
        pos_x = skel.tree(1).posInd(1);
        pos_y = skel.tree(1).posInd(2);
        pos_z = skel.tree(1).posInd(3);
        rot_x = skel.tree(1).expmapInd(1);
        rot_y = skel.tree(1).expmapInd(2);
        rot_z = skel.tree(1).expmapInd(3);
        
        % "For CMU, the absolute position is saved in a different dimension 
        % than it originally was in; we need to write it back to this 
        % dimension"
        dataBase(:,pos_y) = dataBase(:,6);
        
    case 'mit'      % MIT DATA
        pos_x = skel.tree(1).offset(1);
        pos_z = skel.tree(1).offset(2);
        pos_y = skel.tree(1).offset(3);
        rot_x = skel.tree(1).or(1);
        rot_z = skel.tree(1).or(2);
        rot_y = skel.tree(1).or(3);
    otherwise
        error('Unknown skeleton type');
end


% frame 1
% "Assuming we start all at zero (could put in initial conditions here)"
dataBase(1,rot_y) = 0;    % "could be, for example,  newrep{1}(1,3)"
dataBase(1,pos_x) = 0;
dataBase(1,pos_z) = 0;

  
% "integrate" velocities back up into positions--including "constant" term
yaw = cumsum(vertrotdelta);     yaw = [0; yaw(1:end-1)];
m = sqrt(sum(grounddelta.*grounddelta,2));
al = yaw - atan2(grounddelta(:,2),grounddelta(:,1));
x = cumsum(m.*cos(al));         x = [0; x(1:end-1)];
z = cumsum(m.*sin(al));         z = [0; z(1:end-1)];

% store
dataBase(:,[pos_x,pos_z]) = [x, z];

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function dataBase = BC2expCoords(skel,dataBase,yaw)
% convert from body-centred orientations to exponential maps [??]


% extract relative angles
phi = dataBase(:,1);                    % pitch, I think
theta = dataBase(:,2);                  % roll, I think
psi = mod(yaw,2*pi);                    % yaw \in [0,2*pi]


% components of u
switch skel.type
    case 'acclaim', uy = -cos(theta);
    case 'mit',     uy = -cos(phi);         % "easy"
    otherwise,      error('Unknown skeleton type');
end
uz_over_ux = abs(tan(psi));
magux = sqrt((1 - uy.^2)./(1 + uz_over_ux.^2));
maguz = magux.*uz_over_ux;
signz = 2*(psi < pi) - 1;
signx = 2*((psi < pi/2)|(psi > 3*pi/2)) - 1;
uz = maguz.*signz;
ux = magux.*signx;



% components of v ("means different things for MIT and CMU")
switch skel.type
    case 'acclaim',     vy = -cos(phi);
    case 'mit',         vy = -cos(theta);
    otherwise,          error('Unknown skeleton type');
end

% "now we have a quadratic equation to solve for vz"
a = (ux.^2 + uz.^2);
b = 2*uy.*uz.*vy;
c = ((ux.^2).*(vy.^2) - ux.^2 + (uy.^2).*(vy.^2));

% "There are two solutions to the quadratic equation, and one is right
% depending on the quadrant of psi.  Note the other solution corresponds to
% the vector coming out of the other hip."
vz = (-b + signx.*sqrt(b.^2 - 4*a.*c))./(2*a);



% "The rotation matrix is built up a little differently for MIT and CMU"
switch skel.type
    case 'acclaim'  % CMU DATA
        vx = (uy.*cos(phi) - uz.*vz)./ux;
        U = [ux uy uz];     % out right hip
        V = [vx vy vz];     % straight ahead
        WV = -cross(U,V);
        % "NOTE that the r__ subscripts refer to the order, and not to our
        % definition that we have used in the rest of this code.  Yes -- it
        % is confusing!  Here, x- 1st dimension, y- 2nd dimension, z- 3rd 
        % dimension, regardless of what these dimensions represent."
        Rx = U*eye(3);      % what's the point of this??
        Ry = WV*eye(3);
        Rz = V*eye(3);
        RR = cat(1,permute(Rx,[3,2,1]),permute(Ry,[3,2,1]),permute(Rz,[3,2,1]));
        % "Now for display purposes, we want the exponential map 
        % representation.  So do the conversion, based on the rotation 
        % matrix"
        dataBase(:,skel.tree(1).or) = SO32expCoords(RR)';
        
        
    case 'mit'      % MIT DATA
        vx = (uy.*cos(theta) - uz.*vz)./ux;
        U = [ux uz uy];
        V = [vx vz vy];
        WV = cross(U,V);
        
        % "NOTE that the r__ subscripts refer to the order, and not to our 
        % definition that we have used in the rest of this code.  Yes -- it
        % is confusing!  Here, x- 1st dimension, y- 2nd dimension, z- 3rd 
        % dimension, regardless of what these dimensions represent
        Rx = U*eye(3);  %%%% pointless--eliminate
        Ry = V*eye(3);
        Rz = WV*eye(3);
        RR = cat(1,permute(Rx,[3,2,1]),permute(Ry,[3,2,1]),permute(Rz,[3,2,1]));

        % "Now for display purposes, we want the exponential map 
        % representation.  So do the conversion, based on the rotation 
        % matrix"
        dataBase(:,skel.tree(1).or) = -SO32expCoords(RR)';
        % "NOT SURE why rotmat2expmap returns the negative, but the 
        % negative here seems to work.
        
    otherwise
        error('Unknown skeleton type');
end


end
%-------------------------------------------------------------------------%