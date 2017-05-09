function OoR = outofrange(vars,varmin,varmax,str,VERBOSE) %,mods)
%-------------------------------------------------------------------------%
% Revised: 07/02/14
%   -vectorized inputs (vars)
% Revised: 01/07/13
%   -changed strings 'pos' and 'th' to 'Hand-Position' and 'Joint-Angle'
% Revised: 09/04/12
%   -changed initial if... to check for params.NS=prop rather than
%   params.Ndims=2
% Revised: ??/??/??
%   -simplified the calculations for vis by just transforming them into
%   prop and then recursing
% Revised: 02/14/11
%   -fixed calculations for vis, which are really much more complicated
%   than simply checking the min and max.
% Created: 02/03/11
%   -by JGM
%-------------------------------------------------------------------------%


OoR = 0;


% if strcmp(str,'Hand-Position') && strcmp(params.NS,'Hand-Position')
%
%     % for 2D vars, it's easier to check in prop space!
%     Th = IK2link(vars,params,0);
%     jointInd = strcmp(mods,'Joint-Angle');
%     OoR = outofrange(Th,smin(:,jointInd),smax(:,jointInd),'Joint-Angle',VERBOSE,mods);
%
%     if any(OoR) && VERBOSE
%         fprintf(' (vis was the offender)\n');
%     end
%     m
% else
for iDim = 1:size(vars,2)
    excesses = vars(:,iDim) - varmax(iDim);
    iExcessive = (excesses > 0);
    if any(iExcessive)
        errstr = ['warning: ',str,num2str(iDim),' > ','max by ',...
            num2str(excesses(iExcessive)),'\n'];
        OoR = 1;
    end
    excesses = varmin(iDim) - vars(:,iDim);
    iExcessive = (excesses > 0);
    if any(iExcessive)
        errstr = ['warning: ',str,num2str(iDim),' < ','min by ',...
            num2str(excesses(iExcessive)),'\n'];
        OoR = 1;
    end
    
    if OoR && VERBOSE
        fprintf(errstr);
    end
end

% end


end


%%%%%%%%%%%%%%%%%%%%%
% this hand codes in things like sqrt(2)/2 which really ought to be
% calculated from thmin and thmax.  You might want to fix that.
%%%%%%%%%%%%%%%%%%%%%
%     L1 = params.armlengths(1);;
%     L2 = params.armlengths(2);;
%
%     if var(1)^2 + var(2)^2 > (L1 + L2)^2
%         err = sqrt(var(1).^2 + var(2).^2) - (L1 + L2);
%         OoR = 1;
%     elseif (var(1) < -sqrt(2)/2*(L1 + L2)) &&...
%             ((var(1) + sqrt(2)/2*L1)^2 + (var(2) - sqrt(2)/2*L1)^2 > L2^2);
%         err = sqrt((var(1)+sqrt(2)/2*L1)^2 + (var(2)-sqrt(2)/2*L1)^2) - L2;
%         OoR = 1;
%     elseif (var(1) > -sqrt(2)/2*(L1 + L2)) &&...
%             ((var(1)-sqrt(2)/2*L1)^2 + (var(2)-sqrt(2)/2*L1)^2 <= L2^2);
%         err = L2 - sqrt((var(1)-sqrt(2)/2*L1)^2 + (var(2)-sqrt(2)/2*L1)^2);
%         OoR = 1;
%     elseif (var(1) < sqrt(2)/2*(L1 - L2)) &&...
%             (var(1)^2 + var(2)^2 < (L1^2 + L2^2));
%         err = sqrt(L1^2 + L2^2) - sqrt(var(1)^2 + var(2)^2);
%         OoR = 1;
%     end