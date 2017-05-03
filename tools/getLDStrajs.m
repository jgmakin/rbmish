function obsvs = getLDStrajs(params)
% This function basically exists b/c the data may appear either as an array
% (usually from real neural data, where the trajectories are not guaranteed
% to be the same length) or as a structure containing tensors
%
% You may also want to use it to generate quick and dirty data, rather than
% creating a new case in setParams.m and going through generateData.m.

%-------------------------------------------------------------------------%
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%


switch params.datatype
    case 'default'
        
        % stimuli
        %%%%%
        % FIX ME!
        %%%%%
        [S,filterdata] = trajectoryGen(params.dynamics.T,params);
        
        % Ns
        [Ncases,Nstates] = size(filterdata(1).states);
        T = params.dynamics.T;
        Nobsvs = params.N;
        
        % states
        Zcell = mat2cell(shiftdim(cat(3,filterdata(:).states),1),...
            Nstates,T,ones(1,Ncases));
        [obsvs(1:Ncases).Z] = Zcell{:};
        
        % observations
        Y = sampleT(exp(S(:,:)),params.typeUnits{1},params.numsUnits,params);
        Ycell = mat2cell(shiftdim(shortdata(Ncases,3,Y),1),Nobsvs,T,ones(1,Ncases));
        [obsvs(1:Ncases).Y] = Ycell{:};
            
            
    otherwise
        
        
        % generate fresh data
        %%%% broken
        LDSdata = getLDSdata(params);
        %%%%
        
        % if data are in tensors, then put them into an array of mats
        if ndims(LDSdata(1).Z)==3 && length(LDSdata)==1
            
            [Ncases,Nobsvs,T] = size(LDSdata.Y);
            Nstates = size(LDSdata.Z,2);
            
            % observations
            Ycell = mat2cell(shiftdim(LDSdata.Y,1),Nobsvs,T,ones(1,Ncases));
            [obsvs(1:Ncases).Y] = Ycell{:};
            
            % states
            Zcell = mat2cell(shiftdim(LDSdata.Z,1),Nstates,T,ones(1,Ncases));
            [obsvs(1:Ncases).Z] = Zcell{:};
            
            % observation covariance
            if isfield(LDSdata,'SigmaYX')
                SigmaYcell = mat2cell(shiftdim(LDSdata.SigmaYX,1),...
                    Nobsvs,Nobsvs,T,ones(1,Ncases));
                [obsvs(1:Ncases).SigmaYX] = SigmaYcell{:};
            end
            
            % controls
            if isfield(LDSdata,'U')
                Ucell = mat2cell(shiftdim(LDSdata.U,1),Nstates,T,ones(1,Ncases));
                [obsvs(1:Ncases).U] = Ucell{:};
            end
        else
            obsvs = LDSdata;
        end
        
end



end