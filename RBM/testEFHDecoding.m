function [yvar,tErravg] = testEFHDecoding(vishid,hidbiases,visbiases,...
    yvar,testData,params)
% testeEFHDecoding      Test EFH by decoding from it
%
% USAGE:
%   [yvar,tErravg] = testEFHDecoding(vishid,hidbiases,visbiases,...
%        testData,params);
%
% Self-explanatory.  Called by EFH.m
%
% NB: This is not currently rigged up to work on DBNs.  That would require
% a little more thought---tho' not much

%-------------------------------------------------------------------------%
% Cribbed: 01/15/15
%   from EFH.m
%   by JGM
%-------------------------------------------------------------------------%


% prepare the wts
wts{1} = [vishid; hidbiases];
wts{2} = [vishid'; visbiases'];
%%% restarts = ...

% test differently depending on whether it's an EFH or rEFH
if isfield(params,'dynamics')
    [~,~,pEFH] = EFHfilter(testData,wts,params);
    pEFH.name = 'rEFH';
    EFHstats = testDynamics(testData,params,0,pEFH);
    Cvrn = EFHstats(strcmp(params.mods,params.NS)).Cvrn;
else
    NSind = strcmp(params.mods,params.NS);
    [~,R1] = updown(testData.R,wts,params,'means');
    [Shat, ~] = GTPNsuffstats(R1,params);
    Cvrn = cov(Shat(:,:,NSind)-testData.S(:,:,NSind));
end
yvar = [yvar; det(Cvrn)];
tErravg = det(Cvrn);

ydata = gather(yvar);
save('results/EFHtest.mat','ydata');

end