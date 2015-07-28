function paramDisplay(params)

fprintf('\n\nPretraining a %i-layer DBN...\n',length(params.numsUnits));
fprintf('   ')
for i = 1:length(params.numsUnits)-1
    fprintf('%i (%s) x ',params.numsUnits(i),params.typeUnits{i});
end
fprintf('%i (%s)...\n',params.numsUnits(end),params.typeUnits{end});
fprintf('   with %i step(s) of contrastive divergence...\n',...
    params.numCDsteps);
fprintf('   for at most %i epochs...\n',params.DBNmaxepoch);
fprintf('\n\n');    

end

