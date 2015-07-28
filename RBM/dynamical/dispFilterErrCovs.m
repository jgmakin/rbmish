function stats = dispFilterErrCovs(t0,LDSdata,params,varargin)
% THIS FUNCTION APPEARS TO BE RETIRED
% 

%-------------------------------------------------------------------------%
% Revised: 05/06/14
%   -changed the inputs
%   -added line to check for rEFH and use CZ rather than S as the benchmark
%   since CZ doesn't wrap (this is all for the torus case)
% Revised: 04/16/14
%   -cleaned up the tagging
%   -changed dispErrCovs to work better with this, made the associated
%       changes here
% Revised: 01/06/13
%   -now prints figures for all "modalities"
% Revised: 12/16/13
%   -X -> S and etc.
% Created: ??/??/??
%   -by JGM
%-------------------------------------------------------------------------%

% init
S = LDSdata.S;
[Ntraj,Ndims,Nmods,T] = size(S);
samples = t0:T;
q = length(varargin);




% loop through modalities (probably PROP and CTRL)
for iMod = 1:Nmods
    figure(86); clf; hold on;
    
    
    % use unwrapped stimuli for case 'wrapping'
    if strcmp(params.dynamics.walls,'wrapping')
        switch params.mods{iMod}
            case 'Joint-Angle'
                M = params.dynamics.C;
                Ltnts = LDSdata.Z;
            case 'Efference-Copy'
                M = params.dynamics.H;
                Ltnts = LDSdata.U; %%% unfinished! see getLDSdata
        end
        Snotwrapped = zeros(size(S));
        for iTraj = 1:Ntraj
            Snotwrapped(iTraj,:,:) = M*squeeze(Ltnts(iTraj,:,:));
        end
        S(:,:,iMod,:) = Snotwrapped;
    end
    

    
    
    for fltr = 1:q
        % error statistics
        Shat = varargin{fltr}.Shat;
        err = longdata(Shat(:,:,iMod,samples) - S(:,:,iMod,samples));
        stats{fltr}{iMod}.mu = mean(err);
        stats{fltr}{iMod}.cov = cov(err);
        
        % tags for dispErrCovs
        name = varargin{fltr}.name;
        stats{fltr}{iMod}.tags.name = name;
        switch name
            case {'Hand-Position','Joint-Angle','Efference-Copy'}
                stats{fltr}{iMod}.tags.src = 'single';
                stats{fltr}{iMod}.tags.epist = 'empirical';
            otherwise
                stats{fltr}{iMod}.tags.src = 'multiple';
                stats{fltr}{iMod}.tags.epist = varargin{fltr}.name;
        end
        stats{fltr}{iMod}.tags.mod = params.mods{iMod};
        stats{fltr}{iMod}.tags.space = 'neutral';
        stats{fltr}{iMod}.tags.var = 'error';
        
        
        % look at the spatial distribution of errors
        if 0
            switch Ndims
                case 1
                    figure(86); hold on;
                    foo = squeeze(S(:,:,iMod,samples));
                    clr = getColor(stats{fltr}{iMod}.tags.name);
                    scatter(foo(:),err,[],clr);
                    pause();
                case 2
                    figure(86); hold on;
                    foo = longdata(S(:,:,iMod,samples));
                    clr = getColor(stats{fltr}{iMod}.tags.name);
                    % quiver(foo(:,1),foo(:,2),err(:,1),err(:,2),'color',clr)
                    scatter(err(:,1),err(:,2),[],clr);
                    pause();
                otherwise
                    error('unexpected number of stimulus dimensions -- jgm');
            end
        end
        
    end
end


% don't bother plotting if there's only one filter
if q > 1, dispErrCovs(stats,Ntraj*length(samples),params); end

end











