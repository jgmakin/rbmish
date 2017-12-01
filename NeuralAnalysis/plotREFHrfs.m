function plotREFHrfs(wts,params)
% plotREFHrfs   For the rEFH_spikecounts model
%
% USAGE:
%   load([getdir('data'),'RBMish/BMI/wts_rEFH_spikecounts_160407_M1S1_CD1.mat'],'wts','params');
%   plotREFHrfs(wts,params)
%

%-------------------------------------------------------------------------%
% Revised: 09/12/16
%   -functionized
% Created: 07/??/16
%   by JGM
%-------------------------------------------------------------------------%


% fix the params and then get data
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
[X,Q] = params.getLatents([],dataclass,[],'all',params,...
    'sequencelength','singlesequence');
[R,Q] = params.getData(X,Q);

% Ns
NXbins = 25;
NYbins = 25;
% Nunitsperside = 21; % 15;
Nunitsperside = 15;

% vars
[iX,iY,iXdot,iYdot,iXddot,iYddot] = getKinematicStateInds(params.datafile);    
Xddot   = X(:,iXddot);
Yddot   = X(:,iYddot);
Xdot    = X(:,iXdot);
Ydot    = X(:,iYdot);
Y       = X(:,iY);
X       = X(:,iX);
phi     = atan2(Ydot,Xdot);
r       = sqrt(Xdot.^2 + Ydot.^2);

keyboard



% get hidden-unit activities
[~,Z] = updownRDBN(R,wts,params,Q.T);


%%% NB: we make sure that *Y* is in the *rows*
[PrZis1givenYX,MI_ZwithYX] = getCondProbs(Y,X,Z,NYbins,NXbins);
PrZis1givenYX = reshape(PrZis1givenYX,[NXbins*NYbins,size(PrZis1givenYX,3)]);
figure(102);
showrfs(PrZis1givenYX(:,1:(Nunitsperside^2)));
title('position')
matlab2tikzWrapper([params.datatype,'_',num2str(params.Ncdsteps),'_pos'],figure(102))

[PrZis1givenYX,MI_ZwithYX] = getCondProbs(Ydot,Xdot,Z,NYbins,NXbins);
PrZis1givenYX = reshape(PrZis1givenYX,[NXbins*NYbins,size(PrZis1givenYX,3)]);
figure(103);
showrfs(PrZis1givenYX(:,1:(Nunitsperside^2)));
title('velocity')
matlab2tikzWrapper([params.datatype,'_',num2str(params.Ncdsteps),'_vel'],figure(103))

[PrZis1givenYX,MI_ZwithYX] = getCondProbs(Yddot,Xddot,Z,NYbins,NXbins);
PrZis1givenYX = reshape(PrZis1givenYX,[NXbins*NYbins,size(PrZis1givenYX,3)]);
figure(104);
showrfs(PrZis1givenYX(:,1:(Nunitsperside^2)));
title('acceleration')
matlab2tikzWrapper([params.datatype,'_',num2str(params.Ncdsteps),'_acc'],figure(104))

[PrZis1givenYX,MI_ZwithYX] = getCondProbs(phi,r,Z,NXbins,NYbins);
PrZis1givenYX = reshape(PrZis1givenYX,[NXbins*NYbins,size(PrZis1givenYX,3)]);
figure(105);
showrfs(PrZis1givenYX(:,1:(Nunitsperside^2)));
title('speed vs. movement angle')
matlab2tikzWrapper([params.datatype,'_',num2str(params.Ncdsteps),'_spdmvtang'],figure(105))


phiEdges = getEdgeVector(phi,50,[],'3STD');
[PrAis1givenB,MI_AwithB] = getPrAis1givenB(Z,phi,phiEdges);

keyboard


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function [iX,iY,iXdot,iYdot,iXddot,iYddot] = getKinematicStateInds(datafile)

switch datafile
    case {...
            'Indy_datafiles/spikes_and_kinematics_20160407_02.mat',...
            'Indy_datafiles/spikes_and_kinematics_20160411_01.mat',...
            'Indy_datafiles/spikes_and_kinematics_20160627_01.mat',...
            }
        iX = 1; iY = 2; iXdot = 3; iYdot = 4; iXddot = 5; iYddot = 6;
        
        
    case 'Indy_datafiles/spikes_and_kinematics_02_20160407.mat'
        iX = 2; iY = 3; iXdot = 5; iYdot = 6; iXddot = 8; iYddot = 9;
    case 'Indy_datafiles/spikes_and_kinematics_01_20160411.mat'
        iX = 2; iY = 3; iXdot = 8; iYdot = 9; iXddot = 14; iYddot = 15;
    otherwise
        error('bad datafile request -- jgm');
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [PrZis1givenXY,MI_ZwithXY] = getCondProbs(X,Y,Z,NXbins,NYbins)

XEdges = getEdgeVector(X,NXbins,[],'3STD');
YEdges = getEdgeVector(Y,NYbins,[],'3STD');

[PrZis1givenXY,MI_ZwithXY] =...
    getPrAis1givenBandC(Z,X(:),Y(:),XEdges,YEdges);
%[PrAis1givenB,MI_AwithB] = getPrAis1givenB(longdata(Z0),X(:),5*XEdges);

end
%-------------------------------------------------------------------------%