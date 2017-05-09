function J = littleJacobian(mJ,gst)
% littleJacobian    such and such Jacobian from manipulator Jacobian

%-------------------------------------------------------------------------%
% Created: 09/07/10
%   by JGM
%-------------------------------------------------------------------------%

%%%
% vectorized me when you get a chance
%%%

for iCol = 1:size(mJ,2)
    dergstderth = wedge(mJ(:,iCol))*gst;
    J(:,iCol) = dergstderth(1:3,end);
end
    