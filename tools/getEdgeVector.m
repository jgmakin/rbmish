function edgevector = getEdgeVector(X,Nxbins,params,limittype)
% getEdgeVector     Get histogram-bin edges from data 
%
% USAGE:
%   edgevector = getEdgeVector(X,Nxbins,params,'3STD');
%
% getEdgeVector calculates edges for Nxbins bins from data X according to 
% the method specified by the string limittype.  It's especially useful for
% empirical probability or mutual-information calculations (with, e.g.,
% getPrAis1givenB.m or getPrAis1givenBandC.m).
%
% Unless the limittype 'PARAMS' is chosen, the third argument, params, can
% be left empty.  When it is used, it is in standard JGM-EFH format (see
% e.g. setParams.m).  
%
% '3STD' places the penultimate bins at the three-standard-deviation marks,
% and then the final edges at +/- infinity.  Good for normally distributed
% data.
% 
% 'MINMAX' places the penultimate bins just after and before (resp.) the
% min and max of X, and then the final edges at +/- infinity.  Good for
% uniformly distributed data.

%-------------------------------------------------------------------------%
% Cribbed: 07/29/16
%   from getOptimalLag.m
%   by JGM
%-------------------------------------------------------------------------%

fprintf('\n\nHey!! Getting limits with method %s!!\n\n\n',limittype);

switch limittype
    case 'PARAMS'
        N = params.N;
        NSmin = params.smin(1,strcmp(params.NS,params.mods));
        NSmax = params.smax(1,strcmp(params.NS,params.mods));
        Xmin = NSmin;
        Xmax = N/(N-1)*(NSmax - NSmin) + NSmin;
        edgevector = linspace(Xmin,Xmax,Nxbins+1);
        
    case '3STD'
        Xmin = -3*sqrt(var(X(:)));
        Xmax = 3*sqrt(var(X(:)));
        if min(X(:)) > Xmin, Xmin = min(X(:)) + eps; end
        if max(X(:)) < Xmax, Xmax = max(X(:)) - eps; end
        edgevector = linspace(Xmin,Xmax,Nxbins+1-2); % -2 b/c of infs below
        edgevector = [-inf edgevector inf];
        
    case 'MINMAX'
        Xmin = min(X(:)) + eps;
        Xmax = max(X(:)) - eps;
        edgevector = linspace(Xmin,Xmax,Nxbins+1-2); % -2 b/c of infs below
        edgevector = [-inf edgevector inf];
    otherwise
        error('unrecognized limittype -- jgm');
end

end