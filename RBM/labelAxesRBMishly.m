function labelAxesRBMishly(modStr)

switch modStr
    case 'Joint-Angle'
        xlabel('$\prop_1$ (shoulder angle, rad)','Interpreter','none');
        ylabel('$\prop_2$ (elbow angle, rad)','Interpreter','none');
    case 'Efference-Copy'
        fprintf('YOU NEVER FILLED OUT THIS CASE! -- JGM');
    case 'Gains'
        xlabel('$\gain_1$','Interpreter','none');
        ylabel('$\gain_2$','Interpreter','none');
    otherwise
        fprintf('no axis labels! -- jgm\n');
end

end
