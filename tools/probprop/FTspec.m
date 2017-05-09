% HMMspec.m
%
%-------------------------------------------------------------------------%
% Revised: 03/09/12
%   -stuff
% Created: 03/07/08
%   by JGM
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%%
clear all; clc;

% init params
T = 9;
A = [.8 .2; .3 .7];
B = [.8 .2; .4 .6];
p0 = [1; 0];
observations = 'TTFTFFFTT';

% set up some numbering scheme (it hardly matters, as long as no number is
% used twice)
F1 = 1:2:2*T;
X  = 2:2:2*T;
F2 = 2*T + (1:2:2*T);
Y  = 2*T + (2:2:2*T);

% build HMM
tree = cell(T,1);
for i = 1:T
    
    tree{X(i)}.alphabet = {'q','r'};         % hidden
    tree{Y(i)}.alphabet = {'T','F'};         % observed
    
    switch i
        case 1
            tree{F1(i)}.nbhd = [X(i)];
            tree{F1(i)}.fctr = {logical(p0)};% this chooses 1 as the root
            %%% using 'logical' makes sure sumprodFT does the right thing..
            tree{X(i)}.nbhd = [F1(i),F1(i+1),F2(i)];
        case T
            tree{F1(i)}.nbhd = [X(i-1),X(i)];
            tree{F1(i)}.fctr = {A,A'};
            tree{X(i)}.nbhd = [F1(i),F2(i)];
        otherwise
            tree{F1(i)}.nbhd = [X(i-1),X(i)];
            tree{F1(i)}.fctr = {A,A'};
            tree{X(i)}.nbhd = [F1(i),F1(i+1),F2(i)];
    end
    tree{F2(i)}.nbhd = [X(i),Y(i)];
    tree{F2(i)}.fctr = {B,B'};
    tree{Y(i)}.nbhd = [F2(i)]; % ,F3(i)];

    tree{X(i)}.obs = '-';                      % no observation
    tree{Y(i)}.obs = observations(i);
    
end

set(0,'RecursionLimit',2*(T+10))

treeSUM = sumprodFT(tree);


% print the results
for i = [X Y]
    [y ind] = max(exp(treeSUM{i}.logprobs));
    fprintf('%4.2f %6.2f\n',exp(treeSUM{i}.logprobs));
    foo(i) = treeSUM{i}.alphabet{ind};
end
foo
observations
%-------------------------------------------------------------------------%




%-------------------------------------------------------------------------%
%%
clear;


T = 9;
A = [.2 .6 .2; .2 .4 .4; .3 .3 .4];
B = [.3 .7; .4 .6; .3 .7];
p0 = [1; 0; 0];
observations = 'TTFTFFFTT';
% for i = 1:T
%     switch rand > 0.5
%         case 0
%             observations(i) = 'T';
%         case 1
%             observations(i) = 'F';
%     end
% end


% set up some numbering scheme (it hardly matters, as long as no number is
% used twice)
F1 = 1:2:2*T;
X  = 2:2:2*T;
F2 = 2*T + (1:2:2*T);
Y  = 2*T + (2:2:2*T);
set(0,'RecursionLimit',T+10)

% build HMM
tree = cell(T,1);
for i = 1:T
    
    tree{X(i)}.alphabet = {'q','r','s'};     % hidden
    tree{Y(i)}.alphabet = {'T','F'};         % observed
    
    switch i
        case 1
            tree{F1(i)}.nbhd = [X(i)];
            tree{F1(i)}.fctr = {logical(p0)};  % this chooses 1 as the root
            tree{X(i)}.nbhd = [F1(i),F1(i+1),F2(i)];
        case T
            tree{F1(i)}.nbhd = [X(i-1),X(i)];
            tree{F1(i)}.fctr = {A,A'};
            tree{X(i)}.nbhd = [F1(i),F2(i)];
        otherwise
            tree{F1(i)}.nbhd = [X(i-1),X(i)];
            tree{F1(i)}.fctr = {A,A'};
            tree{X(i)}.nbhd = [F1(i),F1(i+1),F2(i)];
    end
    tree{F2(i)}.nbhd = [X(i),Y(i)];
    tree{F2(i)}.fctr = {B,B'};
    tree{Y(i)}.nbhd = [F2(i)]; % ,F3(i)];
    
    tree{X(i)}.obs = '-';                      % no observation
    tree{Y(i)}.obs = observations(i);
    
    % kind of hacky to put these in the first node, but where else??
    tree{1}.fctrinds = [F1 F2];
    tree{1}.rvinds = [X Y];
end

set(0,'RecursionLimit',2*(T+10))


% run sum-product algorithm
treeSUM = sumprodFT(tree);

% print the results
strHID = '%4.2f %6.2f %6.2f\n';
strVIS = '%4.2f %6.2f\n';
for i = [X Y]
    [y ind] = max(exp(treeSUM{i}.logprobs));
    if sum(X==i); str = strHID; else str = strVIS; end
    fprintf(str,exp(treeSUM{i}.logprobs));
    foo(i) = treeSUM{i}.alphabet{ind};
end


foo
observations
