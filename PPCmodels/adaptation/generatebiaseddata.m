function [D0,S0,shft] = generatebiaseddata(Nbatches,shft,params)
% See recalibrationCorePP, the calling function

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -x0 -> S0, with dimension changes
% Created: ??/??/??
%   -by JGM
%-------------------------------------------------------------------------%


% get shft
% shft = getshft(nBatches,stddevs,params);
fprintf('...imposing an initial disrepancy of ');
fprintf('[%01.02f;%01.02f] radians\n',shft(1),shft(2));

% pick a point in the middle
p0.mu = scalefxn([0.5 0.5],[0;0],[1;1],params.roboparams.thmin,params.roboparams.thmax);
p0.cov = 0;

% generate data
[D0, S0] = generateData(Nbatches*params.Ncases,params,'propbias',shft,...
    'stimulusprior',p0);

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function shft = getshft(nBatches,stddevs,params)
%%% this function has the misfortune of generating rather different shifts
%%% for different input gains.  That means that the ICs for the adaptations
%%% under different gains are themselves different, whereas you would
%%% really like them to be the same.


[D0,S0] = generateData(nBatches,params);
SINSMerr = covInCalc(D0,S0,params);
ErrCovIn = SINSMerr{1}.cov + SINSMerr{2}.cov;
iDirection = 1; nDirections = 8;
phi = 2*pi*iDirection/nDirections;
R = [cos(phi) -sin(phi); sin(phi) cos(phi)];
shft = sqrtm(ErrCovIn)*R*[stddevs; 0];

end
%-------------------------------------------------------------------------%
