function [BPTTgradsigW,BPTTgradsigbh] =...
    BPTT(pvisstates,phidmeans,negderLtderzbartr,Wuz,HIDFXN)
% Backpropagation through time (BPTT)

%-------------------------------------------------------------------------%
% Revised: 11/09/15
%   -replaced derLtderbhtr with its negative, changed signs accordingly
% Revised: 07/20/15
%   -replaced two arguments with their difference(s), derLtderbhtr.
% Created: 05/29/15 (happy b'day, VMO)
%   by JGM
%-------------------------------------------------------------------------%

% Ns
[T,Nhid] = size(phidmeans);

% the derivative of the feedforward function
if length(HIDFXN) > 1
    error('you never programmed in this case - jgm');
else
    switch HIDFXN{1}
        case 'Bernoulli'
            fprimevec = @(zhat,tt)(zhat(tt,:).*(1-zhat(tt,:)));
        otherwise
            error('you haven''t implemented these cases yet -- jgm\n');
    end
end
    
% malloc (matrix whose *rows* are g_t' [see paper])
G = zeros(T,Nhid,'like',pvisstates);
% Assume the final prediction is perfect, since you don't actually have
% the sensory information at time T+1.  Hence the zeros.
G(end,:) = zeros(1,Nhid,'like',pvisstates).*fprimevec(phidmeans,T);
for t = (T-1):-1:1
    G(t,:) = (negderLtderzbartr(t+1,:) + G(t+1,:)*Wuz).*fprimevec(phidmeans,t);
end
BPTTgradsigW = pvisstates'*G;
BPTTgradsigbh = sum(G,1);

end


