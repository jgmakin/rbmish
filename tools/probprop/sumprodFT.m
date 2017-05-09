function tree = sumprodFT(tree)
% SUMPRODFT() sum-product algorithm for factor trees
%
% Hand this function a factor tree (see FTspec.m for an example tree
% specification) and it'll return the marginal probabilities.  The 
% algorithm follows M.I. Jordan's pseudocode in Chapter 4 of his
% (still-unpublished) _An Introduction to Probabilistic Graphical Models_.
%
% Don't try this on really long chains (~4000 if there are three hidden
% states and two observed states), or the recursion will crash matlab.
%
% NB. THERE'S A MINOR ISSUE WITH ZERO PROBABILITIES: B/c you work with log
% probabilies to avoid precision errors, p=0 gives logp=-Inf, which causes
% a problem since (-Inf - -Inf) ~= 0, as you might think, but to NaN.  To
% obviate this, issue, you've set evidence and the prior (which happens to
% be, essentially, evidence, being [1; 0]) to *logicals*; the code then
% looks for such and treats them differently, sc. not putting them into log
% space until after the problematic operation.

%-------------------------------------------------------------------------%
% Revised: 03/13/12
%   -completed: yields the same results as HMMspec on the equivalent FT
% Revised: 03/02/12
%   -fixed all sorts of problems
% Created: 02/22/12
%   by JGM
%-------------------------------------------------------------------------%

% init
f = 1;                                      % set root node
%%% this should prolly be set by the HMM...

tree = Evidence(tree);                      % apply observations
tree = muCollect(0,f,tree);                 % recursively collect
tree = nuDistribute(0,f,tree);              % recursively distribute
tree = ComputeMarginalL(tree);              % compute marginal probabilties


end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function tree = Evidence(tree)

for i = 1:length(tree)
    if ~isfield(tree{i},'fctr')
        delta = strcmp(tree{i}.obs,tree{i}.alphabet);
        if sum(delta);                          % if there be evidence
            N = length(tree);
            tree{N+1}.nbhd = i;
            tree{N+1}.fctr = {delta(:)};        % set the evidence potential
            % does this say that there's only evidence at the leaves???
            tree{i}.nbhd = [tree{i}.nbhd N+1];
        end
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = nuCollect(s,i,tree)

% fprintf('%i is a factor node\n',s);
nbhd = tree{i}.nbhd;
nbhdnots = nbhd(nbhd~=s);
for t = nbhdnots
    tree = muCollect(i,t,tree);
end
tree = nuSendMessageL(i,s,tree);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = muCollect(i,s,tree)

% fprintf('%i is a rv node\n',i);
nbhd = tree{s}.nbhd;
nbhdnoti = nbhd(nbhd~=i);
for j = nbhdnoti
    tree = nuCollect(s,j,tree);
end
tree = muSendMessageL(s,i,tree);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = nuSendMessageL(i,s,tree)
% send a message *from* i *to* s
% I'm a RANDOM-VARIABLE node!
% NB: messages are log(messages).

% note that all the neighbors are factor nodes
nbhd = tree{i}.nbhd;
nbhdnots = nbhd(nbhd~=s);
m = 0;
for t = nbhdnots
    tnbhd = tree{t}.nbhd;
    m = m + tree{t}.mu{tnbhd==i};
    % the message from k to j
end

% get sum
if s > 0
    tree{i}.nu{nbhd==s} = m;
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = muSendMessageL(s,i,tree)
% send a message *from* s *to* i
% I'm a FACTOR node!
% NB: messages are log(messages).

MAXEXPONENT = log(realmax);

% note that all the neighbors are RV nodes
nbhd = tree{s}.nbhd;
nbhdnoti = nbhd(nbhd~=i);
m = 0;
for j = nbhdnoti
    jnbhd = tree{j}.nbhd;
    m = m + tree{j}.nu{jnbhd==s};
    % the (nu) message from j to s
end

% get sum
if i>0
    fctr = tree{s}.fctr{nbhd==i};
    if islogical(fctr)
        b = log(double(fctr)) + repmat(m',size(fctr,1),1);
        tree{s}.mu{nbhd==i} = sum(b,2);  % shouldn't do anything
    else
        b = log(fctr) + repmat(m',size(fctr,1),1);
        B = MAXEXPONENT - log(length(m)) - max(b,[],2);
        tree{s}.mu{nbhd==i} = log(sum(exp(b+repmat(B,1,size(b,2))),2)) - B;
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = nuDistribute(i,s,tree)

if i > 0
    tree = nuSendMessageL(i,s,tree);
end
nbhd = tree{s}.nbhd;
nbhdnoti = nbhd(nbhd~=i);
for j = nbhdnoti
    tree = muDistribute(s,j,tree);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = muDistribute(s,i,tree)

tree = muSendMessageL(s,i,tree);
nbhd = tree{i}.nbhd;
nbhdnots = nbhd(nbhd~=s);
for t = nbhdnots
    tree = nuDistribute(i,t,tree);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function tree = ComputeMarginalL(tree)

MAXEXPONENT = log(realmax);
    
for i = 1:length(tree)
    
    if ~isfield(tree{i},'fctr')
        nbhd = tree{i}.nbhd;

        % "multiply" the outgoing and incoming messages
        s = nbhd(1);                          % any edge will do!!
        m = tree{s}.mu{tree{s}.nbhd==i} + tree{i}.nu{nbhd==s};
        
        % normalize these unnormalized log probabilities
        b = m;
        B = MAXEXPONENT - log(length(m)) - max(b);
        normalizer = log(sum(exp(b+B))) - B;
        tree{i}.logprobs = b - normalizer;
    
    end
    
end

end
%-------------------------------------------------------------------------%

