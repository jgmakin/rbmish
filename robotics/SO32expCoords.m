function omegath = SO32expCoords(R)
% SO32expCoords  Convert from SO(3) to "exponential coordinates"
%
% USAGE:
%   omegath = SO32expCoords(R);
% We sometimes (e.g., the mocap data set) need "exponential coordinates" of
% a rotation, i.e. omega*th, with omega the axis, and th the amount, of the
% rotation.  We'd like to convert to this product from a rotation matrix.
%
% N.B.: The map from SO(3) to exponential coordinates is *many-to-one*.
% Even when restricting the range of th \in [0,2*pi), there are two
% distinct solutions for (omega,th)--and for R = I, an infinite number of
% solutions (any axis will do, since there is no rotation).  This function
% simply chooses based on the output of matlab's acos; and in the case of
% R=I, sets the rotation axis to [1;0;0];
%
% See Murray, Li, Sastry pp. 29-30 for details.
% 
% The input must have size 3 x Nexamples (but Nexamples may be 1) 
%
% See also the inverse function, expCoords2SO3.m. 

%-------------------------------------------------------------------------%
% Created: 11/22/16
%   by JGM
%-------------------------------------------------------------------------%


% Ns
Nexamples = size(R,3);

% th = acos( (Tr(R)-1)/2 )
%%% until matlab allows arrayfun on gpuarrays...
%%% TrR = arrayfun(@(ii)(trace(R(:,:,ii))),1:Nexamples);
%%% th = acos( (TrR-1)/2 );
TrR = arrayfun(@(ii)(trace(R(:,:,ii))),1:Nexamples,'UniformOutput',0);
th = acos( ([TrR{:}]-1)/2 );
iNoRot = abs(th) < eps;
%%%

% omega = [r32-r32; r13-r31; r21 - r12]/(2*sin(th))
num = permute(cat(1,...
    R(3,2,:) - R(2,3,:),...
    R(1,3,:) - R(3,1,:),...
    R(2,1,:) - R(1,2,:)),[1,3,2]);
num(:,iNoRot) = repmat([1; 0; 0],[1,sum(iNoRot)]);
den = 2*sin(th);
den(:,iNoRot) = 1;
omega = num./den;
    
% multiply (w.l.o.g because omega must have length 1)
omegath = omega.*th;

end



























