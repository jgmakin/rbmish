function ressum = checkKFSteadyState(CvrnMat,params)
% This function compares the predicted steady-state covariance matrix
% (after a *time update*) to the predicted (ditto) covariance matrix, by
% solving the appropriate discrete-time "Riccati equation."  It also checks
% to make sure the system is controllable and observable, b/c otherwise the
% Riccati equation is not solvable.

%-------------------------------------------------------------------------%
% Created: 04/30/13
%   by JGM
%-------------------------------------------------------------------------%


% params
A = params.A;
C = params.C;
SigmaX = params.SigmaX;
SigmaYX = params.SigmaYX;
Nstates = size(A,1);


% solve the discrete-time algebraic Riccati equation for the steady
% state covariance
if rank(ctrb(A,SigmaX)) < Nstates
    fprintf('uncontrollable...\n');
    ressum = NaN;
elseif rank(obsv(A,C)) < Nstates
    fprintf('unobservable...\n');
    ressum = NaN;
else
    [P,L,G] = dare(A',C',SigmaX,SigmaYX);
    % are P and the final time-updated covariance equal?
    res = P - (A*CvrnMat(:,:,end)*A' + SigmaX);
    ressum = sum(abs(res(:)));
end

% plot the covariance about velocity
% figure(12); clf; plot(1:size(CvrnMat,3),squeeze(CvrnMat(2,2,:)))


end
%-------------------------------------------------------------------------%