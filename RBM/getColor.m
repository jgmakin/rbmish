function clr = getColor(name)

setColors
switch name
    case 'Hand-Position'
        clr = VIScolor;
    case 'Joint-Angle'
        clr = PROPcolor;
    case 'Gaze-Angle'
        clr = EYEcolor;
    case 'Efference-Copy'
        clr = ECcolor;
    case {'rEFH','EFH'}
        clr = EFHcolor;
    case 'opt'
        clr = OPTcolor;
    case 'EM'
        clr = EMcolor;
    case 'EM (best)'
        clr = EMBESTcolor;
    case 'obs'
        clr = OBScolor;
    case {'obs (no u)','obs 1st-ord'}
        clr = OBSNOUcolor;
    case 'prior'
        clr = [1 0 0];
    otherwise
        fprintf('warning: unknown estimator name -- jgm\n');        
        clr = [1 0 0];
end


% if strcmp(tags.src,'single')                % input
%     switch tags.mod
%         case 'Hand-Position'
%             clr = VIScolor;
%         case 'Joint-Angle'
%             clr = PROPcolor;
%         case 'Gaze-Angle'
%             clr = EYEcolor;
%         case 'Efference-Copy'
%             clr = ECcolor;
%         otherwise
%             error('unknown modality -- jgm');
%     end
% else                                        % integ
%     switch tags.epist
%         case {'empirical','rEFH','EFH'}
%             clr = EFHcolor;
%         case {'theoretical','opt'}
%             clr = OPTcolor;
%         case 'EM'
%             clr = EMcolor;
%         case 'EM (best)'
%             clr = EMBESTcolor;
%         case 'obs'
%             clr = OBScolor;
%         case {'obs (no u)','obs 1st-ord'}
%             clr = OBSNOUcolor;
%         otherwise
%             error('unrecognized source of knowledge -- jgm');
%     end
% end

end