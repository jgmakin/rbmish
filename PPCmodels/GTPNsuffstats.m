function [cntrsOfMass, ttlSpks] = GTPNsuffstats(R,params)
% GTPNsuffstats     Sufficient statistics of Gaussian-tuned Poisson neurons
%
%   USAGE:
%       [cntrsOfMass, ttlSpks] = GTPNsuffstats(R,params);
%
% Given a matrix of GTPN populations (Nexamples x Nmods*N^Ndims) and the
% params structure, compute the centers of mass (Nexamples x Ndims x Nmods)
% and the total spike counts (Nexamples x Nmods).

%-------------------------------------------------------------------------%
% Revised: 07/09/14
%   -transposed output ttlSpks
% Revised: 07/07/14
%   -rewrote from scratch to work on tensors rather than using loops
% Revised: 01/03/14
%   -fixed addition from last revision
% Revised: 12/18/13
%   -now works on data PPCs containing more than one "modality."  NB: this
%   changes the shape of the structure fields, Shat and CvrnMat!!!
% Revised: 07/02/13
%   -added input argument tuningCov, which improves the logic of the
%   function (it shouldn't have itself to decide to use the prop cov).
% Renamed: 07/02/13
%   -from PPC.m to GTPNsuffstats.m
% Cribbed: 06/28/13
%   from KF4PPC
%   by JGM
%-------------------------------------------------------------------------%

% Ns
mods = params.mods;
Nmods = length(mods);
Ndims = params.Ndims;
N = params.N;
Nexamples = size(R,1);
Nlattice = N^Ndims;

% (Nexamples x Nmods*N^Ndims) -> (Nmods*Nexamples x N^Ndims)
Rtall = reshape(permute(reshape(R,[Nexamples,Nlattice,Nmods]),[3,1,2]),...
    [Nmods*Nexamples,Nlattice]);
ttlSpksTall = Rtall*ones(Nlattice,1,'like',R);

% construct "preferred directions" on an Ndims lattice
latticePDs = ndlattice(Ndims,ones(1,1,'like',R):N);

% for data on the Ndims-torus, center the data first
gridShifts = zeros(Nmods*Nexamples,Ndims);
if strcmp(params.walls,'wrapping')
    [Rtall,gridShifts] = centerNoisyHills(Rtall,params.N);
end

% calculate center of mass
cntrsOfMassTall = (Rtall*latticePDs)./ttlSpksTall - gridShifts;

% put back into standard tensor shapes (Nexamples x Ndims x Nmods)
cntrsOfMass = permute(reshape(cntrsOfMassTall,[Nmods,Nexamples,Ndims]),[2 3 1]);
ttlSpks = reshape(ttlSpksTall,[Nmods,Nexamples])';

% in case of no spikes!
for iMod = 1:Nmods
    theseCoM = cntrsOfMass(:,:,iMod);
    badCoMInds = isnan(theseCoM(:,1));  %%% only need to check dim=1;
	badCoMInds = gather(badCoMInds);   
 
    if any(badCoMInds)
        fprintf('warning: found %d bad centers of mass ',sum(badCoMInds));
        fprintf('in %s; replacing with guesses -- jgm\n',mods{iMod});
        
        % replace with...other, good centers of mass!
        goodCoM = theseCoM(~badCoMInds,:);
        randInds = ceil(rand(sum(badCoMInds),1)*sum(~badCoMInds));
        theseCoM(badCoMInds,:) = goodCoM(randInds,:);
        cntrsOfMass(:,:,iMod) = theseCoM;
    end
end

%%%%%% might consider adding this
% fprintf('\n\ndoing terrible clamping-at-edges thing to avoid complex numbers\n\n');
% smin = shiftdim(params.smin,-1);
% smax = shiftdim(params.smax,-1);
% cntrOfMass = (cntrOfMass >= smin).*cntrOfMass + (cntrOfMass < smin).*smin;
% cntrOfMass = (cntrOfMass <= smax).*cntrOfMass + (cntrOfMass > smax).*smax;
%%%%%%


% convert into real-world units
for iMod = 1:Nmods
    cntrsOfMass(:,:,iMod) = grid2world(cntrsOfMass(:,:,iMod),...
        [params.smin(:,iMod) params.smax(:,iMod)],params);
end



end
