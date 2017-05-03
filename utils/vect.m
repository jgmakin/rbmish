function v = vect(T)
% vect  Vectorize multi-dimensional objects (matrices, tensors, etc.)
%
% v = vect(T) just performs v = T(:). The reason for this function is that,
% in practice, T is often created by some function, and this avoids having
% to create a temporary variable in the current workspace.

%-------------------------------------------------------------------------%
% Created: 10/06/16
%   by JGM
%-------------------------------------------------------------------------%

v = T(:);

end