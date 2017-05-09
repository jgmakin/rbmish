function u = unifSmplAboutCenter(cntrs,swing,Nsamples)
% unifSmplAboutCenter   Sample uniformly about a center point
%
% USAGE
%   u = unifSmplAboutCenter(cntrs,swing,Nsamples)
%   g = unifSmplAboutCenter(gains,params.swing,Nexamples)
%
% Given a vector of centers (1 x Ncntrs), a "swing" (fraction of cntr that 
% the random variable can swing above or below it), and the requested 
% number of samples, produces a matrix (Nsamples x Ncntrs) of samples whose
% columns are each uniformly distributed cntr(iCol) with minimum
% (cntr - swing*cntr) and maximum (cntr + swing*cntr).  

%-------------------------------------------------------------------------%
% Created: 06/21/14
%   by JGM
%-------------------------------------------------------------------------%

Ncntrs = size(cntrs,2);
u = (1 + swing*2*(rand(Nsamples,Ncntrs,'like',cntrs)-0.5)).*cntrs;

end