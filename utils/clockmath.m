function t3 = clockmath(binaryop,arg1,arg2)
% clockmath     Perform addition and subtraction on times
%
% USAGE:
%   t3 = clockmath(@minus,[0,10,35,29],[0,09,46,37])
%   t3 = clockmath(@plus,[0,0,48,52],[0,09,46,37; 0,11,24,10])
%   t3 = clockmath(@times,[0,0,7,1],4.62)
%   t3 = clockmath(@rdivide,[0,0,32.0000,24.4000],4.68);
%
% It's very annoying to have to keep converting seconds to minutes and so
% on, so this function takes care of that clock math.  At the moment it's
% extremely rudimentary but could be expanded to cover more cases.  
%
% NB: The input (and output) times must be 4-vectors containing:
%
%   [day,hour,minute,second]
%
% You could imagine expanding this to dates.

%-------------------------------------------------------------------------%
% Revised: 10/07/16
%   -modified to allow matrix time arguments (i.e., multiple times at once)
% Created: 10/05/16
%   by JGM
%-------------------------------------------------------------------------%

% bases of the different "digits"
basevec = [365,24,60,60];

% the naive operation
t3a = bsxfun(binaryop,arg1,arg2);

% for borrowing and carrying
t3b = t3a;
carryborrowvec = floor(t3b./basevec);
while any(carryborrowvec(:)~=0)
    t3b = mod(t3b,basevec) +...
        [carryborrowvec(:,2:end),zeros(size(carryborrowvec,1),1)];
    %%% at the moment, changes in *years* are ignored
    carryborrowvec = floor(t3b./basevec);
end

% for division and multiplication: we want integer units everywhere---
%   except the final "digit" (seconds), where decimals can't be helped
t3c = floor(t3b);
fractionalpart = (t3b - t3c).*basevec;
t3d = t3c + [zeros(size(fractionalpart,1),1),fractionalpart(:,1:end-1)];

% output
t3 = t3d;


end