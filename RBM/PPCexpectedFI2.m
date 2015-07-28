function FI = PPCexpectedFI2(range,g,params,xtrue)
% PPCexpectedFI     Expected Fisher information, Pois or Bino population
%   PPCEXPECTEDFI computes the Fisher information matrix FI for a
%   population of Gaussian-tuned neurons, encoding a parameter in RANGE and
%   with additional parameters specified by PARAMS.
%-------------------------------------------------------------------------%
% Created: 12/14/10
%   by JGM
%-------------------------------------------------------------------------%


% init params
Ndims = params.Ndims;                               % encoded vars/neuron
N = params.N;
C = params.C;
n = params.nexperiments;
respLength = params.respLength;
gridsize0 = params.gridsize;
r = range(:,2) - range(:,1);
[Ncases,twom,Nbatches] = size(xtrue);

% rescale (if necessary)
tuningCov = (diag(r(:)/respLength)*sqrtm(C))^2;
gridsize = scalefxn(gridsize0*ones(Ndims,1),zeros(Ndims,1),...
    respLength*ones(Ndims,1),zeros(Ndims,1),r);
S = inv(tuningCov);


% average over the actual x
% init
J = 0;
b = zeros(N,Ndims);
for i = 1:Ndims
    b(:,i) = linspace(0,gridsize(i),N)' + range(i,1);
end
tic
% loop through stimuli (v)
for iCase = 1:Ncases
    for iBatch = 1:Nbatches
        
        % compute f
        v = xtrue(iCase,:,iBatch);
        if Ndims == 1
            d = S*(v - b).^2;
            f = exp(-d/2);
            Ja = g*f.*(S^2*(b - v).^2);
            if strcmp(params.typeUnits{1},'Binomial')
                Ja = Ja.*(n./(n-g*f));
            end
            % Ja = g*(n./(n-g*f)).*f.*(S^2*(b - v).^2);
            J = J + sum(Ja(:));
        elseif Ndims == 2      

            z1z1 = repmat((v(1) - b(:,1)).^2,1,N);
            z1z2 = (v(1) - b(:,1))*(v(2) - b(:,2)');
            z2z2 = repmat((v(2) - b(:,2)').^2,N,1);
            
            s1z1z1 = S(1,1)*z1z1;
            s0z1z2 = S(1,2)*z1z2;
            s4z2z2 = S(2,2)*z2z2;
            s0s0z1z1 = S(1,2)*S(2,1)*z1z1;
            s1s4z1z2 = S(1,1)*S(2,2)*z1z2;
            s0s0z2z2 = S(1,2)*S(2,1)*z2z2;
            
            s1s1z1z1 = S(1,1)*s1z1z1;
            s1s0z1z1 = S(1,2)*s1z1z1;
            s1s0z1z2 = S(1,1)*s0z1z2;
            s0s0z1z2 = S(1,2)*s0z1z2;
            s4s0z1z2 = S(2,2)*s0z1z2;
            s4s0z2z2 = S(1,2)*s4z2z2;
            s4s4z2z2 = S(2,2)*s4z2z2;
        
            d = s1z1z1 + 2*s0z1z2 + s4z2z2;
            f = exp(-d/2);
            Ja11 = f.*(s1s1z1z1 + 2*s1s0z1z2 + s0s0z2z2);
            Ja12 = f.*(s1s0z1z1 + s1s4z1z2 + s0s0z1z2 + s4s0z2z2);
            Ja22 = f.*(s0s0z1z1 + 2*s4s0z1z2 + s4s4z2z2);
            if strcmp(params.typeUnits{1},'Binomial')
               Ja11 = (n./(n-g*f)).*Ja11;
               Ja12 = (n./(n-g*f)).*Ja12;
               Ja22 = (n./(n-g*f)).*Ja22;
            end
            Ja(1,1) = g*sum(sum(Ja11));
            Ja(1,2) = g*sum(sum(Ja12));
            Ja(2,1) = Ja(1,2);
            Ja(2,2) = g*sum(sum(Ja22));
            
            J = J + Ja;
        else
            error('This fxn not written for Ndims > 2\n -- jgm');
        end
    end
end
toc
FI = J/Ncases/Nbatches;

end