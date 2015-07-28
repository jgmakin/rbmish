function tree = maxprod(tree)
% MAXPROD() max-product algorithm for trees
%
% Hand this function a tree (see HMMspec.m for an example tree
% specification) and it'll return the most likely configuration the vars as
% well as the probability of that configuration.  The algorithm follows 
% M.I. Jordan's pseudocode in Chapter 4 of his (still-unpublished) 
% _An Introduction to Probabilistic Grapsical Models_.


%-------------------------------------------------------------------------%
% Revised: 03/02/12
%   -changed a lot of things; most imp., implemented numerically safer
%   computations by using logprobs instead of probs.
% Created: 02/22/12
%   by JGM
%-------------------------------------------------------------------------%

% init
f = 1;                                      % set root node
%%% this should prolly be set by the HMM...

tree = Evidence(tree);                      % apply observations
tree = Collect(0,f,tree);
tree = Distribute(0,f,tree);


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function tree = Evidence(tree)

for i = 1:length(tree)
    delta = strcmp(tree{i}.obs,tree{i}.alphabet);
    if sum(delta);                          % if there be evidence
        tree{i}.psiS = delta(:);            % set the evidence potential
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = Collect(i,j,tree)

nbhd = tree{j}.nbhd;
nbhdnoti = nbhd(nbhd~=i);
for k = nbhdnoti
    tree = Collect(j,k,tree);
end
tree = SendMessageL(j,i,tree);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = SendMessageL(j,i,tree)
% send a message *from* j *to* i
% NB msgs are now log(messages).

nbhd = tree{j}.nbhd;
nbhdnoti = nbhd(nbhd~=i);
m = zeros(size(tree{j}.alphabet))';
for k = nbhdnoti
    knbhd = tree{k}.nbhd;
    m = m + tree{k}.msgs{knbhd==j};
    % the message from k to j
end

% get max
psiS = tree{j}.psiS;    
foo = (m + log(double(psiS)));
if i > 0
    psi = tree{j}.psi{nbhd==i};
    [Y I] = max(log(psi) + repmat(foo',size(psi,1),1),[],2);
else
    [Y I] = max(foo);
end

% assign to this node
tree{j}.msgs{nbhd==i} = Y;                  % max
tree{j}.cnfg = I;                           % argmax

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = Distribute(i,j,tree)

tree = SetValue(i,j,tree);      % tree.xstar...
nbhd = tree{j}.nbhd;
nbhdnoti = nbhd(nbhd~=i);
for k = nbhdnoti
    tree = Distribute(j,k,tree);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = SetValue(i,j,tree)

if i > 0
    cnfg = tree{j}.cnfg;
    tree{j}.xstar = cnfg(tree{i}.xstar);
else
    tree{j}.xstar = tree{j}.cnfg;
end

end
%-------------------------------------------------------------------------%






%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function tree = SendMessage(j,i,tree)
% send a message *from* j *to* i
% NB msgs are now log(messages).


nbhd = tree{j}.nbhd;
nbhdnoti = nbhd(nbhd~=i);
m = 1;
for k = nbhdnoti
    knbhd = tree{k}.nbhd;
    m = m.*tree{k}.msgs{knbhd==j};
    % the message from k to j
end

% get max
psiS = tree{j}.psiS;
if i > 0
    psi = tree{j}.psi{nbhd==i};
    [Y I] = max(psi*diag(psiS.*m),[],2);
else
    [Y I] = max(psiS.*m);
end

% assign to this node
tree{j}.msgs{nbhd==i} = Y;                  % max
tree{j}.cnfg = I;                           % argmax

end
%-------------------------------------------------------------------------%