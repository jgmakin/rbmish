function [Rout,Z] = updownRDBN(Rin,wts,params,T,varargin)
% updownRDBN    updownDBN for recurrent DBNs
%
% USAGES:
%       
%   Rout = updownRDBN(Rin,wts,params,T)
%
%   [Rout,Z] = updownRDBN(Rin,wts,params,T)
% 


%-------------------------------------------------------------------------%
% Created: 01/05/17
%   by JGM
%-------------------------------------------------------------------------%

%%% TO DO:
% (1) allow samples as well as means


% params
numsUnits = params.numsUnits;
typeUnits = params.typeUnits;

% Ns
Ntraj   = size(Rin,1)/T;
hidNums = numsUnits{end};
Nlayers = length(wts)/2+1;          % in *this*, possibly partial, DBN
Mlayers = length(params.numsUnits); % in the *complete* DBN

% init/malloc
Rin     = shortdata(Ntraj,3,Rin);
Rout    = zeros(size(Rin),'like',Rin);
Z       = zeros([Ntraj,sum(hidNums),T],'like',Rin);
thisZ   = zeros([Ntraj,sum(hidNums)],'like',Rin);
VERBOSE = 1;
propagation = 'means';
for i=1:length(varargin)
    switch varargin{i}
        case {'suffstats','MLE','modes','means','Nsamples'}
            propagation = varargin{i};
        case 'quiet'
            VERBOSE = 0;
        otherwise
            fprintf('unrecognized option for updownRDBN -- jgm\n');
    end
end
if VERBOSE
    fprintf('propagating ');
    switch propagation
        case 'Nsamples'
            fprintf('%i sample(s)',params.smpls);
        case 'means'
            fprintf('means');
        otherwise
            error('failed: unrecognized method for updownRDBN!');
    end
    fprintf(' up through DBN\n');
end

% up-down pass
fprintf('\nFiltering (EFH)');
for t = 1:T
    
    % push data up to the penultimate layer
    thisR = Rin(:,:,t);
    for iLayer = 1:length(wts)
        
        % current layer info
        jLayer  = Nlayers-abs(Nlayers-iLayer-1); % really!
        dstrbs  = typeUnits{jLayer};
        nums    = numsUnits{jLayer};
        
        % if you're in the penultimate layer, cat recurrent activities
        if jLayer == Mlayers, thisR = [thisZ,thisR]; end
        
        % if you're in the ultimate layer, output only the right half
        if iLayer == Mlayers
            thisZ = thisR;
            thisR = invParamMap(thisR,...
                wts{iLayer}(1:end-1,(end-sum(nums)+1):end),...
                wts{iLayer}(end,(end-sum(nums)+1):end),dstrbs,nums,params);
        else
            thisR = invParamMap(thisR,wts{iLayer}(1:end-1,:),...
                wts{iLayer}(end,:),dstrbs,nums,params);
        end
        
        % sample?
        switch propagation
            case 'means'
                % do nothing
            case 'Nsamples'
                if iLayer<jLayer % only sample on the way up
                    NNNN = params.smpls;
                    states = sampleT(thisR,dstrbs,nums,params);
                    for i=1:NNNN-1
                        states = states + sampleT(thisR,dstrbs,nums,params);
                    end
                    thisR = states/NNNN;
                end
        end
        
        
    end
    Z(:,:,t)    = thisZ;
    Rout(:,:,t) = thisR;
    
end
fprintf('\n');

Rout = longdata(Rout);
Z = longdata(Z);


end