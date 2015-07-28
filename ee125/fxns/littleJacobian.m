function J = littleJacobian(mJ,gst)
% LITTLEJACOBIAN    such and such Jacobian from manipulator Jacobian
%-------------------------------------------------------------------------%
% Created: 09/07/10
%   by JGM
%-------------------------------------------------------------------------%

for iCol = 1:size(mJ,2)
    dergstderth = wedge(mJ(:,iCol))*gst;
    J(:,iCol) = dergstderth(1:3,end);
end
    