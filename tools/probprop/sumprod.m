function tree = sumprod(tree)
% SUMPROD() sum-product algorithm for trees
%
% Hand this function a tree (see HMMspec.m for an example tree
% specification) and it'll return the marginal probabilities.  The 
% algorithm follows M.I. Jordan's pseudocode in Chapter 4 of his
% (still-unpublished) _An Introduction to Probabilistic Graphical Models_.
%
% Don't try this on really long chains (~4000 if there are three hidden
% states and two observed states), or the recursion will crash matlab.

%-------------------------------------------------------------------------%
% Revised: 03/02/12
%   -fixed all sorts of problems
% Created: 02/22/12
%   by JGM
%-------------------------------------------------------------------------%

% init
f = 1;                                      % set root node
%%% this should prolly be set by the HMM...

tree = Evidence(tree);                      % apply observations
tree = Collect(0,f,tree);
tree = Distribute(0,f,tree);
tree = ComputeMarginalL(tree);


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
% NB messages are now log(messages).

MAXEXPONENT = log(realmax);

nbhd = tree{j}.nbhd;
nbhdnoti = nbhd(nbhd~=i);
m = 0;
for k = nbhdnoti
    knbhd = tree{k}.nbhd;
    m = m + tree{k}.msgs{knbhd==j};
    % the message from k to j
end

% get sum
if i > 0
    foo = m + log(double(tree{j}.psiS));
    psi = tree{j}.psi{nbhd==i};
    b = log(psi) + repmat(foo',size(psi,1),1);   
    B = MAXEXPONENT - log(length(m)) - max(b,[],2);
    tree{j}.msgs{nbhd==i} = log(sum(exp(b+repmat(B,1,size(b,2))),2)) - B;
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = Distribute(i,j,tree)

if i > 0
    tree = SendMessageL(i,j,tree);
end
nbhd = tree{j}.nbhd;
nbhdnoti = nbhd(nbhd~=i);
for k = nbhdnoti
    tree = Distribute(j,k,tree);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = ComputeMarginalL(tree)

MAXEXPONENT = log(realmax);
    
for i = 1:length(tree)
    
    nbhd = tree{i}.nbhd;
    m = 0;
    for j = nbhd
        jnbhd = tree{j}.nbhd;
        m = m + tree{j}.msgs{jnbhd==i};
        % the message from k to j
    end
    
    b = log(double(tree{i}.psiS)) + m;
    B = MAXEXPONENT - log(length(m)) - max(b);
    normalizer = log(sum(exp(b+B))) - B;
    tree{i}.logprobs = b - normalizer;
    
end

end
%-------------------------------------------------------------------------%






%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
function tree = ComputeMarginal(tree)


for i = 1:length(tree)
    
    nbhd = tree{i}.nbhd;
    m = 1;
    for j = nbhd
        jnbhd = tree{j}.nbhd;
        m = m.*tree{j}.msgs{jnbhd==i};
        % the message from k to j
    end
    
    probs = (tree{i}.psiS).*m;
    tree{i}.probs = probs/sum(probs);
end

end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function tree = SendMessage(j,i,tree)
% send a message *from* j *to* i

nbhd = tree{j}.nbhd;
nbhdnoti = nbhd(nbhd~=i);
m = 1;
for k = nbhdnoti
    knbhd = tree{k}.nbhd;
    m = m.*tree{k}.msgs{knbhd==j};
    % the message from k to j
end

% get sum
psiS = tree{j}.psiS;
if i > 0
    psi = tree{j}.psi{nbhd==i};
    tree{j}.msgs{nbhd==i} = sum(psi*diag(psiS.*m),2); % psi*(psiS.*m);  % sum
end

end
%-------------------------------------------------------------------------%