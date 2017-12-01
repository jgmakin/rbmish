function are_equal = tensorEq(T1,T2)
% tensorEq  Check if two tensors are equal
%
% USAGE:
%   sum_abs = tensorEq(T1,T2)
% 
% Often you want to check if two tensors are equal.  It's annoying to type
% all the composed commands, so here it is in one step

%-------------------------------------------------------------------------%
% Created: 09/09/17
%   by JGM
%-------------------------------------------------------------------------%

are_equal = ~sum(abs(vect(T1-T2)));

end