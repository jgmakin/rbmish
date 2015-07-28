function [BPTTgradsigW,BPTTgradsigbh] =...
    BPTT(pvisstates,phidmeans,derLtderbhtr,Wuz,HIDFXN)
% Backpropagation through time (BPTT)

%-------------------------------------------------------------------------%
% Revised: 07/20/15
%   -replaced two arguments with their difference(s), derLtderbhtr.
% Created: 05/29/15 (happy b'day, VMO)
%   by JGM
%-------------------------------------------------------------------------%

% Ns
[T,Nhid] = size(phidmeans);

% the derivative of the feedforward function
switch HIDFXN
    case 'Bernoulli'
        fprimevec = @(zhat,tt)(zhat(tt,:).*(1-zhat(tt,:)));
    otherwise
        error('you haven''t implemented these cases yet -- jgm\n');
end

% malloc (matrix whose *rows* are g_t'*F_t, F_t the jacobian of f)
GF = zeros(T,Nhid,'like',pvisstates);
% Assume the final prediction is perfect, since you don't actually have
% the sensory information at time T+1.  Hence the zeros.
GF(end,:) = zeros(1,Nhid,'like',pvisstates).*fprimevec(phidmeans,T);
for t = (T-1):-1:1
    GF(t,:) = ((derLtderbhtr(t+1,:) + GF(t+1,:))*Wuz).*fprimevec(phidmeans,t);
end
BPTTgradsigW = -pvisstates'*GF;
BPTTgradsigbh = -sum(GF,1);

end

%%% Sutskever's version
% GF(end,:) = zeros(1,Nhid,'like',pvisstates).*fprimevec(phidmeans,T);
% for t = (T-1):-1:1
%     derLPENderbhtr = qhidstates(t,:) - phidstates(t,:);
%     GF(t,:) = (derLPENderbhtr + GF(t+1,:))*Wuz.*fprimevec(phidmeans,t);
% end
% BPTTgradsigW = -pvisstates(1:end-1,:)'*GF(2:end,:);
% BPTTgradsigW((Nhid+1):end,:) = BPTTgradsigW((Nhid+1):end,:)*Wuz;
% BPTTgradsigbh = -sum(GF,1);
%%%