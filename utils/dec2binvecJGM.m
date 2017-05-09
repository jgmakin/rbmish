function bin = dec2binvecJGM(dec,n)
% dec2binvecJGM     Convert from decimal to binary (vector) representation
%
% USAGE:
%   bin = dec2binvecJGM(dec,n);
% 
% Here n is the number of bits; without it, the minimal number of bits is
% used.
%
% There are two reasons to rewrite this function:
%   (1) dec2bin.m, which outputs a *string* rather than a vector, is not in
%   a special toolbox; but (remarkably) dec2binvec.m is!
%   (2) dec2bin.m works on vector inputs (returning matrix outputs), but
%   (remarkably), dec2binvec.m does not!
%
% NB that the *least* significant difit is listed *first*.

%-------------------------------------------------------------------------%
% Cribbed: 12/29/16
%   from matlab's dec2bin.m.
%-------------------------------------------------------------------------%

if nargin < 2, n = 1; end
    
% How many digits do we need to represent the numbers?
[~,e]=log2(max(dec(:)));
bin = fliplr(rem(floor(dec(:)*pow2(1-max(n,e):0)),2));

end


