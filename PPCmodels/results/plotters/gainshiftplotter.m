% compute fisher info in whichever population is "neutral"
% [D0,S0] = generateData(1000,params);
% SINSMerr = covInCalc(D0,S0,params);
% clear D0;
% ErrCovIn = SINSMerr{1}.cov + SINSMerr{2}.cov;


%-------------------------------------------------------------------------%
% Revised: 12/12/13
%   -changed of statsN based on new output for estStatsCorePP.m
% Created: ??/??/??
%   by JGM
%-------------------------------------------------------------------------%


nDirections = 8;
nRadii = 6;
setColors
iShift = 0;
NSIND = find(strcmp(params.mods,params.NS));
Nexamples = 1000;

%
close all
% for iShift = 1:nRadii
    % for iDirection = 0% :nDirections-1
        for i = 1:length(gainratios) 
            
            % i = iDirection+nDirections*(iShift-1)+2;
        
        
            if iDirection == 0
                dispErrCovs(ErrorStatsArray(i,[1,3:4]),40000,params,...
                    2*iShift);
                %%%% HACK: writing in 40000 explicitly here is only
                %%%% "correct" b/c you (think you) used nbatches = 1000 as
                %%%% the argument to generateData inside test.m.  You always
                %%%% do, but this should still be fixed
                
                % axis([biasOPT(1,i)-0.001 biasOPT(1,i)+0.004,...
                %   biasOPT(2,i)-0.001 biasOPT(2,i)+0.011]);
                
                % get this shift
                % shft = shiftArray(:,i);
                
                % generate input and output data
                gains = mean([params.gmin; params.gmax]);
                params.gmin = gains;
                params.gmax = gains;
                params.smpls = 15;
                  
                [Di,Si] = generateData(Nexamples,params,'propbias',shft);              
                Do = updownDBN(Di,wts,params,'Nsamples');
                [statsL,statsN,ShatL,ShatN] = estStatsCorePP(Si,params,'CoM',Do);
                
                
                errors = ShatN(:,:,NSIND) - Si(:,:,NSIND);
                offset = ErrorStatsArray{i,3}{NSIND}.mu - statsN{NSIND}.mu;
                shftErrs = errors + repmat(offset',size(ShatN,1),1);
                
                hold on;
                scatter(shftErrs(:,1),shftErrs(:,2),5,EFHcolor,'x');
                hold off;
            end
            
        
        end
     % end
% end

% for i = 1:2:13
%     figure(i);
%     close
% end