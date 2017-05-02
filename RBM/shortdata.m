function varargout = shortdata(Nrows,Nways,varargin)
% shortdata     Transform from matrix to "standard" EFH-training tensor
%
% USAGES:
%   Nways = 4; Nrows = 30;
%   Sshort = shortdata(Nrows,Nways,Slong);
%
%   Nways = 3; Nrows = 30;
%   [Rshort,Vshort] = shortdata(Nrows,Nways,Slong,Vlong);
% 
% Vshort = shortdata(Nrows,Nways,Vlong) undoes Vlong = longdata(Vshort);
% i.e., it transforms the (Nexamples x Ndims x N3 x N4 x ...) tensor V into
% a "short" (along the first way) tensor with Nrows in the first "way" (the 
% rows) and Nexamples/Nrows in the last one.  So, e.g., when N3=N4=...=1,
% Nways = 3 forces Vshort to have size (Nrows x Ndims x Nexamples/Nrows).
% Nways = 4 produced Vshort w/size (Nrows x Ndims x N3 x Nexamples/Nrows).
% 
% NB that it is reasonable to request Nways=4 (or higher) even when N3=1
% (and so on for higher ways), but NOT conversely; i.e., it is not sensible
% to ask for Nways = 3 when N3 > 1.
% 
% NB! Consecutive examples in Vmat (along its first way) become consecutive
% examples along the *final* way of Vtensor!!
% 
% See also longdata.m.
%
% Canonical usages:
%   Sshort = shortdata(Ncases,4,Slong);
%       (Nexamples x Ndims x Nmods) -> (Ncases x Ndims x Nmods x Nbatches)
%
%   Rshort = shortdata(Ncases,3,Rlong);
%       (Nexamples x Nunits -> (Ncases x Nunits x Nbatches)


%-------------------------------------------------------------------------%
% Revised: 12/07/16
%   -rewrote help, renamed the variables
% Revised: 12/18/13
%   -added an input argument for the "order" of the tensor.  This
%   (unfortunate) argument is needed to cover the case where length(params.
%   Nmods)=1 and therefore the longdata'd input of stimuli S has only order
%   2, and therefore looks like a data matrix (D), whereas its shortform 
%   should have order 4, not 3.  IF YOU CAN THINK OF SOMETHING BETTER, then
%   please change this.
% Revised: 12/10/13
%   -changed to work with tensors with more than three "dimensions"---as
%   long as the first and last are the number of cases and batches!
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%

Nargs = nargin - 2;
for i = 1:Nargs
    longform = varargin{i};
    
    longdims = size(longform);
    cmptdims = [longdims(1)/Nrows Nrows longdims(2:end)];
    
    if rem(longdims(1),Nrows) ~= 0
        error('number of cases (input 1) incompatible with Nexamples');
    else
        
        cmptform = reshape(longform,cmptdims);       
        shrtform = permute(cmptform,[2:Nways,1]);
        varargout{i} = shrtform;

    end
end

end
