function subdirs = getSubdirs
% There's a bug in Matlab2014a (and perhaps earlier) that prevents
% dir('*.') from returning directories with names longer than eight
% characters (presumably only on Windows machines).  This is a workaround.

%-------------------------------------------------------------------------%
% Xmitted: 04/16/14
%   by 
%       Shridhar Shah
%       MathWorks Technical Support Department 
%   to 
%       JGM
%-------------------------------------------------------------------------%

A = dir;
myDirs = [A(:).isdir];
subdirs = {A(myDirs).name}';
subdirs(ismember(subdirs ,{'.','..'})) = [];

end