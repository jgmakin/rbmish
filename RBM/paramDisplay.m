function paramDisplay(params)

%-------------------------------------------------------------------------%
% Revised: 02/24/16
%   -rewrote to accommodate new typeUnits, numsUnits format
% Created: ??/??/??
%   -by JGM
%-------------------------------------------------------------------------%

fprintf('\n\nPretraining a %i-layer DBN...\n\t',length(params.numsUnits));

% the first layer
nums = params.numsUnits{1};
dstrbs = params.typeUnits{1};
fprintf('Nunits: \t\t\t [%i (%s)',nums(1),dstrbs{1});
for iGrp = 2:length(dstrbs)
    fprintf(' + %i (%s)',nums(iGrp),dstrbs{iGrp});
end
fprintf(']');

% all the remaining layers
for iLayer = 2:length(params.numsUnits)
    nums = params.numsUnits{iLayer};
    dstrbs = params.typeUnits{iLayer};
    
    fprintf(' x [%i (%s)',nums(1),dstrbs{1});
    for iGrp = 2:length(dstrbs)
        fprintf(' + %i (%s)',nums(iGrp),dstrbs{iGrp});
    end
    fprintf(']');
end
fprintf('\n');

% the other things to say
fprintf('\tNcdsteps: \t\t\t %i\n',params.Ncdsteps);
fprintf('\tNcases/batch: \t\t %i\n',params.Ncases);
fprintf('\tNbatches/epoch: \t %i\n',params.Nbatches);
fprintf('\tNepochs: \t\t\t %i\n',params.NepochsMax);
fprintf('\n\n');    

end

