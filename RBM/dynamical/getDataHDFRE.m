function [R,Q] = getDataHDFRE(S,Q,charset,visDstrb)
% getDataHDFRE
% 
% USAGE:
%   R = getDataHDFRE(S,Q,charset,params.typeUnits{1}{1})

%-------------------------------------------------------------------------%
% Created: 12/??/16
%   by JGM
%-------------------------------------------------------------------------%

%%% TO DO:
% (1) modify for GPU


% "Ns"
Nexamples   = size(S,1);
T           = Q.T;
Nseq        = floor(Nexamples/T);
Nchars      = length(charset);

% load the data and remove useless characters
completetext = fileread([getdir('data'),'HDFRE/HDFRE.txt']);
completetext(completetext==sprintf('\n')) = ' ';
ISGOODCHAR = any(completetext == charset');
completetext = completetext(ISGOODCHAR);

% find the beginnings of sentences
capitalletters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
periodinds = find(completetext == '.');
capitalletterinds = find(any(completetext == capitalletters'));
Nsentences = length(periodinds);
[~,bin] = histc(periodinds,capitalletterinds);
sentenceinitinds = capitalletterinds(bin+1);

% now pull out Nseq sequences that start at the beginnings of sentences
randinds = ceil(Nsentences*rand(Nseq,1));
sentenceinitinds = sentenceinitinds(...
    sentenceinitinds < length(completetext)-T);
startinds = sentenceinitinds(randinds);
allinds = startinds' + (0:T-1);
examplesentences = permute(completetext(allinds),[1,3,2]);
Ronehot = longdata(examplesentences) == charset(2:end);

% conversion from strings to binary representation
switch visDstrb %%% assume there's only one
    case 'Categorical'
        R = double(Ronehot);
    case 'Bernoulli'
        bintags = dec2binvecJGM((0:(Nchars-1))',7);
        examplecharsinds = Ronehot*(1:(Nchars-1))' + 1;
        R = bintags(examplecharsinds,:);
    otherwise
        error('unexpected type of units for %s -- jgm',...
            params.datatype);
end


end