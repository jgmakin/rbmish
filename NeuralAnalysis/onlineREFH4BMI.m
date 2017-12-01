function [xhatDynamic,znow,xhatTU,CvrnTU] = onlineREFH4BMI(...
    r,zprev,xhatTU,CvrnTU,Bv,LDSparams,wts,params)
% onlineREFH4BMI    online BMI decoding with the rEFH
%
% USAGE:
%   [xhatDynamic,znow,xhatTU,CvrnTU] =
%       onlineREFH4BMI(r,zprev,xhatTU,CvrnTU,Bv,LDSparams,wts,params);
%
% Essentially, this is a combination of testEFHBMI.m and filterDecoder.m,
% except that the LDSparams are not fit internally but passed in (probably
% having been computed by filterDecoder, in fact).

% Inputs:
%   r           vector of current spike counts
%   zprev       previous hidden-unit activities
%   xhatTU      E[X_t|xhat_{t-1},...,xhat_0]        (time update)
%   CvrnTU      Cvrn[X_t|xhat_{t-1},...,xhat_0]     (time update)
%   Bv          "static" rEFH decoder (assumes *all* units are decoded)
%   LDSparams   params for the linear dynamical system on which KF is run
%   wts         rEFH weights
%   params      rEFH params
%
% Outputs:
%   xhatDynamic E[X_t|xhat_t,...,xhat_0]            (measurement update)
%   znow        current hidden-unit activities
%   xhatTU      E[X_{t+1}|xhat_t,...,xhat_0]        (time update)
%   CvrnTU      Cvrn[X_{t+1}|xhat_t,...,xhat_0]     (time update)
%
% The measurement update xhatDynamic should be used to write to the screen;
% the time updates xhatTU and CvrnTU (along with znow) are just for the 
% next iteration.
%
% NB: There is an unfortunate convention conflict between the KF and the
% rEFH.  For the former, data at a single time point are stored in a column 
% vector; for the latter, they are stored in a row vector.  (This is baked
% into the structures wts and LDSparams, so it's not really feasible to
% resolve the conflict in favor of one or the other.)  Hence, xhat* are 
% column vectors, but r and z* are row vectors (and similarly for the
% parameters).



%-------------------------------------------------------------------------%
% Created: 05/06/16
%   -by JGM
%-------------------------------------------------------------------------%

% Ns
numsUnits = params.numsUnits;
Nlayers = length(numsUnits);

% up and down through the rEFH
%%%
[rhat,znow] = updownRDBN(r,wts,params,1,'verbosity',0,'initial recurrents',zprev);
% If this is too slow, it can replaced with a much simpler, less flexible
% version.
%%%
uhat = invParamMap(znow,wts{Nlayers}(1:end-1,1:sum(numsUnits{end})),...
    wts{Nlayers}(end,1:sum(numsUnits{end})),params.typeUnits{end},...
    params.numsUnits{end},params);
v = [rhat,uhat,znow,1];

% now decode
xhatStatic = Bv'*v';
[xhatDynamic,xhatTU,CvrnTU] = KFonestep(xhatStatic,xhatTU,CvrnTU,...
    LDSparams.A,LDSparams.muX,LDSparams.SigmaX,...
    LDSparams.C,LDSparams.muYX,LDSparams.SigmaYX);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [xhatMU,xhatTU,CvrnTU] = KFonestep(y,xhatTU,CvrnTU,...
    A,muX,SigmaX,C,muYX,SigmaYX)

% measurement update
K = CvrnTU*C'/(SigmaYX + C*CvrnTU*C');
innov = (y - muYX) - C*xhatTU;
xhatMU = xhatTU + K*innov;
CvrnMU = CvrnTU - K*C*CvrnTU;

% time update
CvrnTU = A*CvrnMU*A' + SigmaX;
xhatTU = A*xhatMU + muX;

% % measurement update
% CtrInfo = C'/SigmaYX;
% InfoMU = CtrInfo*C + InfoTU;
% xhatMU = InfoMU\(CtrInfo*(y - muYX) + InfoTU*xhatTU);
% 
% % time update
% InfoTU = inv(A/InfoMU*A' + SigmaX);
% xhatTU = A*xhatMU + muX;


end
%-------------------------------------------------------------------------%

