function Q = mocapBase2EFHdata(skel,Motion)
% mocapBase2EFHdata  Convert motion-capture motion coordinates
%
% USAGE:
%   Motion = mocapBase2EFHdata(skel,Motion)
% 
% Converts the original mocap data (either 'acclaim' or 'mit') into the
% form used by Sutskever2009: "body-centered coordinates and...ground-plane
% differences," and thence into EFH-ready data.

%-------------------------------------------------------------------------%
% Revised: 11/21/16
%   -mostly rewrote: tensorized, functionized, rationalized
%   by JGM
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
% This program preprocesses data (stage 1 of 2)
% We need to unwrap the rotation about the vertical axis
% before we take the difference in the next stage of preprocessing
%
% We support two types of skeletons:
%  1) Those built from the CMU database (acclaim)
%     http://mocap.cs.cmu.edu/
%  2) Those built from data from Eugene Hsu (mit)
%     http://people.csail.mit.edu/ehsu/work/sig05stf/
%
% The program assumes that the following variables are set externally:
% n1     -- order of the first layer CRBM
% Motion -- cell array of motion data
% skel   -- skeleton structure
%-------------------------------------------------------------------------%


% for each of the sequences...
for iSeq = 1:length(Motion)
    
    % "body-centred coordinates"
    [pitchroll, yaw] = motion2RPY(skel,Motion{iSeq});
    
    % "ground-plane vectors"
    [velu, velv] = motion2grndVels(skel,Motion{iSeq},yaw);
    
    % now write out
    % "even if it is CMU data, let's overwrite the first 6 dimensions"
    switch skel.type
        case 'acclaim'  % CMU DATA
            % "Absolute vertical position - dimension 2 is the one we don't
            % want to touch"
            Motion{iSeq}(:,6) = Motion{iSeq}(:,2);
        case 'mit'      % MIT DATA
            % "The absolute vertical is already in dimension 6"
    end
    Motion{iSeq}(:,1:2) = pitchroll;          % body-centred orientations
    Motion{iSeq}(:,3) = [diff(yaw); yaw(end)-yaw(end-1)];
    Motion{iSeq}(:,4:5) = [velu velv];        % ground-plane relative vels
    
end

% now convert the body coordinates to EFH-usable data
Q = bodyCoords2EFHdata(skel,Motion);


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [pitchroll, yaw] = motion2RPY(skel,motion)

% get some useful things for this data type
% "find the rotation matrix based on the root Exponential map angles"
switch skel.type
    case 'acclaim'  % CMU DATA
        rootRotInd = skel.tree(1).expmapInd;
        R0 = [0 1 0; 0 0 1; 1 0 0];
    case 'mit'      % MIT DATA
        rootRotInd = skel.tree(1).or;
        R0 = eye(3);
    otherwise
        error('Unknown skeleton type');
end
Nframes = size(motion,1);

% base orientation
R = expCoords2SO3(motion(:,rootRotInd)');

% CMU data set uses right multiplication
if strcmp(skel.type,'acclaim'), R = permute(R,[2,1,3]); end

% rotate the first two "normals"
R1 = tensorOp(R,repmat(R0(:,1:2),[1,1,Nframes]));

% pitch and roll (rotation of x-y plane wrt the vertical)
pitchroll = permute(acos(...
    tensorOp(repmat(-R0(:,3)',[1,1,Nframes]),R1)./sqrt(sum(R1.*R1))),...
    [3,2,1]);

% yaw (rotation within the base x-y plane)
switch skel.type
    case 'acclaim'
        yaw = squeeze(atan2(R1(3,2,:),R1(1,2,:)));
    case 'mit'
        yaw = squeeze(atan2(R1(2,1,:),R1(1,1,:)));
end
yaw = unwrap(yaw);
% "We now do phase unwrapping to recognize the fact that within a
% sequence we may turn more than 2pi radians. The unwrapping maintains
% the continuity of the turn"


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [velu, velv] = motion2grndVels(skel,motion,yaw)

% get some useful things for this data type
switch skel.type
    case 'acclaim'  % CMU DATA
        grdinds = [1,3];                            % -x, z
    case 'mit'      % MIT DATA
        grdinds = [4,5];                            % -x, y
    otherwise
        error('Unknown skeleton type');
end

% horizontal position in terms of the body-centred coordinate system
horizVel = diff(motion(:,grdinds));
horizSpeed = sqrt(sum(horizVel.*horizVel,2));
mvmtAngle = atan2(horizVel(:,2),horizVel(:,1));

% For MIT, the u (cos) component is the component in line with the
% "straight out of body" vector
% For CMU, the u (cos) component is the compoent in line with the
% "straight out of right hip" vector
velu = horizSpeed.*cos(yaw(1:(end-1)) - mvmtAngle);
velv = horizSpeed.*sin(yaw(1:(end-1)) - mvmtAngle);
velu(end+1)=velu(end);   velv(end+1)=velv(end); % hack (in original)

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function data = bodyCoords2EFHdata(skel,Motion)
% bodyCoords2EFHdata   Convert motion-capture body coords into EFH data 

%-------------------------------------------------------------------------%
% Revised: 11/21/16
%   -funtionized
%   -replaced loops, repmats with bsxfun and arrayfun
%   by JGM
% Cribbed: 11/21/16
%   from preprocess2ilya.m
%-------------------------------------------------------------------------%

switch skel.type
    case 'acclaim'  % CMU-style data
        % No offsets, but several dimensions are constant
        indx = [ 1:6 ...
            10:12 13 16:18 19 ...
            25:27 28 31:33 34 ...
            37:39 40:42 43:45 46:48 49:51 52:54 ...
            58:60 61 65 67:69 73:75 ...
            79:81 82 86 88:90 94:96 ];
        % root (special representation)
        % lfemur ltibia lfoot ltoes
        % rfemur rtibia rfoot rtoes
        % lowerback upperback thorax lowerneck upperneck head
        % (lclavicle ignored) lhumerus lradius lwrist lhand
        %   (fingers are constant) lthumb
        % (rclavicle ignored) rhumerus rradius rwrist rhand
        %   (fingers are constant) rthumb
    case 'mit'      % MIT-style data
        indx = [   1:6 7:9 14 19:21 26 31:33 38 43:45 50 55:57 61:63 67:69 ...
            73:75 79:81 85:87 91:93 97:99 103:105 ];
        
        % save the offsets, they will be inserted later
        data.offsets = [  Motion{1}(1,10:12); Motion{1}(1,16:18); ...
            Motion{1}(1,22:24); Motion{1}(1,28:30); Motion{1}(1,34:36); ...
            Motion{1}(1,40:42); Motion{1}(1,46:48); Motion{1}(1,52:54); ...
            Motion{1}(1,58:60); Motion{1}(1,64:66); Motion{1}(1,70:72); ...
            Motion{1}(1,76:78); Motion{1}(1,82:84); Motion{1}(1,88:90); ...
            Motion{1}(1,94:96); Motion{1}(1,100:102); Motion{1}(1,106:108)];
    otherwise
        error('Unknown skeleton type');
end

% combine the data and z-score
alldata = cat(1,Motion{:});
alldata = alldata(:,indx);
data.mu = mean(alldata);
data.sig = std(alldata);
data.R = (alldata - data.mu)./data.sig;

% build an index of sequences
Nseq = length(Motion);
data.seqindex = mat2cell(1:size(alldata,1),1,arrayfun(@(ii)(size(Motion{ii},1)),1:Nseq));

end
%-------------------------------------------------------------------------%