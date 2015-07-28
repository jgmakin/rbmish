function varargout = shortdata(Ncases,order,varargin)
% Undo longdata.m, i.e. split the examples (which must be in the first
% dimension) into cases and batches (which will be in the first and last
% dimensions, respectively, of the output).
%
% USAGE:
%   Sshort = shortdata(numcases,order,Slong)
%   [Sshort,Yshort] = shortdata(numcases,order,Slong,Ylong)


%-------------------------------------------------------------------------%
% Revised: 12/18/13
%   -added an input argument for the "order" of the tensor.  This
%   (unfortunate) argument is needed to cover the case where params.Nmods=1
%   and therefore the longdata'd input of stimuli S has only order 2, and
%   therefore looks like a data matrix (D), whereas its shortform should
%   have order 4, not 3.  IF YOU CAN THINK OF SOMETHING BETTER, than please
%   change this.
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
    cmptdims = [longdims(1)/Ncases Ncases longdims(2:end)];
    
    if rem(longdims(1),Ncases) ~= 0
        error('number of cases (input 1) incompatible with Nexamples');
    else
        
        cmptform = reshape(longform,cmptdims);       
        shrtform = permute(cmptform,[2:order,1]);
        varargout{i} = shrtform;

    end
end

end
