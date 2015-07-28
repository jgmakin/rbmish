function beta = hiddecoderPP(D0,x0,wts,params,varargin)
% HIDDECODER    Linear regression on the hidden units
%   HIDDECODER returns linear-regression coefficients for decoding the
%   output.  You can request as many coeffecient vectors as there are input
%   modalities (params.Nmods).

%-------------------------------------------------------------------------%
% Created: 02/01/11
%   by JGM
%-------------------------------------------------------------------------%

%%%%% if you really want to keep this, it should be vectorized (which would
%%%%% be very easy) (07/09/14, jgm)


% init
DISP = 0;

% variable input arguments
for i=1:2:length(varargin)
    switch varargin{i}
        case 'display'
            DISP = varargin{i+1};
        otherwise
            fprintf('unrecognized option -- jgm\n');
    end
end

% generate (deepest) hidden unit activities
fprintf('generating hidden unit activities...\n');
[V,x] = getHiddenActs(D0,x0,wts,params);

% solve the normal eq'ns
leftinv = inv(V'*V)*V';
beta = leftinv*x;

% plot
if DISP
    for i = 1:size(x0,2)
        figure();
        scatter(V*beta(:,i),x(:,i));
    end
end

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [V,x] = getHiddenActs(D0,x0,wts,params)

% init
[Ncases,Nvis,Nbatches] = size(D0);
Nexamples = Ncases*Nbatches;
Nhid = params.numsUnits(end);

% malloc
hiddenmeans = zeros(Ncases,Nhid,Nbatches);

% loop
tic
HADBEENCLOSED = isempty(gcp('nocreate'));
if HADBEENCLOSED, pool = parpool(4); end
parfor (iBatch = 1:Nbatches, 8)
    probs = D0(:,:,iBatch);
    for iLayer = 1:length(wts)/2
        HIDFXN = params.typeUnits{iLayer+1};
        probs = feedforward(probs,wts{iLayer}(1:end-1,:),...
            wts{iLayer}(end,:),HIDFXN,params);
    end
    % fprintf('training data fraction: %f\n',batch/size(D0,3));
    hiddenmeans(:,:,iBatch) = probs;
end
if HADBEENCLOSED, delete(pool); end
toc

% put into regression shape
Vv = reshape(shiftdim(hiddenmeans,1),params.numsUnits(end),Nexamples)';
V = [Vv ones(Nexamples,1)];
x = reshape(shiftdim(x0,1),size(x0,2),Nexamples)';

end
%-------------------------------------------------------------------------%