% BMIhelper script
% 
% ADDING RESULTS OF A NEW DECODER TO EXISTING RESULTS
%
% Suppose you run BMImaster on a *new* decoder.  This little script joins
% these new results to the existing file, and saves the result

%-------------------------------------------------------------------------%
% Revised: 10/03/17
%   -changed to accommodate cases where the new data are a *replacement*
%   for old data, rather than a new decoder entirely
% Created: 09/09/17
%   bu JGM
%-------------------------------------------------------------------------%


%%%%
% (1) SHOULD REALLY PUT IN SOME ERROR CHECKERS!  E.G., USING NDATATEST.
%%%%

clear; clc;

% user parameters
monkeys = {'Indy','Loco','Chewie'};
newDecoderName = 'kfemdynamic';
file_suffix = 'WF';
sensory_unit_type = 'Poisson';
swept_param = 'binwidths';
TOPLOT = 1;
TOSAVE = 0;


% vector lengths
Nmonkeys = length(monkeys);

% filename parts
data_dir = sprintf('%sRBMish/BMI/',getdir('data'));
sweep_name = [upper(swept_param(1)), swept_param(2:end-1), 'Sweep'];
for iMk = 1:Nmonkeys
    
    % load the backup version of the old complete results (all decoders)
    mk = monkeys{iMk};
    old_load_file = sprintf('%sRsqs_%s_%s_%s.mat',...
        data_dir,sweep_name,sensory_unit_type,mk);
    old_results = load(old_load_file);
    
    % construct the file name of the new data (under new decoder)
    new_load_file = sprintf('%sassemble_me/Rsqs_%s_%s_%s_%s.mat',...
        data_dir,sweep_name,sensory_unit_type,mk,file_suffix);
    new_results = load(new_load_file);
    
    % append the new decoder to the master list
    iDecoderOld = strcmp(old_results.decodernames,newDecoderName);
    iDecoderNew = strcmp(new_results.decodernames,newDecoderName);
    if any(iDecoderOld)
        if TOPLOT
            for iState = 1:size(old_results.Rsqs,2)
                for iSweep = 1:size(old_results.Rsqs,3)
                    figure(7); clf; hold on;
                    keyboard
                    scatter(old_results.Rsqs(:,iState,iSweep,iDecoderOld),...
                        new_results.Rsqs(:,iState,iSweep,iDecoderNew))
                    title(sprintf('state %i, sweep var %i',iState,iSweep))
                    plot([0,1],[0,1],'k');
                    hold off
                    pause()
                end
            end
        end
        old_results.Rsqs(:,:,:,iDecoderOld) = new_results.Rsqs(:,:,:,iDecoderNew);
    else
        Ndecoders_old = length(old_results.decodernames);
        old_results.decodernames{Ndecoders_old+1} = newDecoderName;
        old_results.Rsqs = cat(4, old_results.Rsqs, new_results.Rsqs);
    end
    
    
    % construct the file name of the composite results
    out_file = sprintf('%sRsqs_%s_%s_%s.mat',...
        data_dir,sweep_name,sensory_unit_type,mk);

    if TOSAVE
        save(out_file, '-struct', 'old_results');
    end
end

