% HMMspec.m
%
%-------------------------------------------------------------------------%
% Created: 3/7/08
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

% build HMM
tree = cell(T,1);
for i = 1:T
    tree{i}.alphabet = {'q','r'};           % hidden
    tree{i+T}.alphabet = {'T','F'};         % observed
    switch i
        case 1
            tree{i}.nbhd = [2,i+T];
            tree{i}.psiS = p0;              % this chooses 1 as the root
            tree{i}.psi = {A',B'};
        case T
            tree{i}.nbhd = [i-1,i+T];
            tree{i}.psiS = [1; 1];          % non-root singletons get 1
            tree{i}.psi = {A,B'};            % pair potentials
        otherwise
            tree{i}.nbhd = [i-1,i+1,i+T];
            tree{i}.psiS = [1; 1];          % non-root singletons get 1
            tree{i}.psi = {A,A',B'};          % pair potentials
    end
    tree{i+T}.nbhd = [i];
    tree{i+T}.psiS = [1; 1];                % non-root singletons get 1
    tree{i+T}.psi = {B};                    % pair potentials
    
    tree{i}.obs = '-';                      % no observation
    tree{i+T}.obs = observations(i);        % 
    
end

% run max-product algorithm
treeMAX = maxprod(tree);

% print the results
for i = 1:T
    configstar(i) = treeMAX{i}.alphabet{treeMAX{i}.xstar};
end
configstar
MAP = treeMAX{1}.msgs{1}



% run max-product algorithm
treeSUM = sumprod(tree);

% print the results
for i = 1:length(treeSUM)
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


set(0,'RecursionLimit',T+10)

% build HMM
tree = cell(T,1);
for i = 1:T
    tree{i}.alphabet = {'q','r','s'};       % hidden
    tree{i+T}.alphabet = {'T','F'};         % observed
    switch i
        case 1
            tree{i}.nbhd = [2,i+T];
            tree{i}.psiS = p0;              % this chooses 1 as the root
            tree{i}.psi = {A',B'};          % pair potentials
        case T
            tree{i}.nbhd = [i-1,i+T];
            tree{i}.psiS = [1; 1; 1];       % non-root singletons get 1
            tree{i}.psi = {A,B'};           % pair potentials
        otherwise
            tree{i}.nbhd = [i-1,i+1,i+T];
            tree{i}.psiS = [1; 1; 1];       % non-root singletons get 1
            tree{i}.psi = {A,A',B'};        % pair potentials
    end
    tree{i+T}.nbhd = [i];
    tree{i+T}.psiS = [1; 1];                % non-root singletons get 1
    tree{i+T}.psi = {B};                    % pair potentials
    
    tree{i}.obs = '-';                      % no observation
    tree{i+T}.obs = observations(i);        % 
    
end


% run max-product algorithm
treeMAX = maxprod(tree);

% print the results
for i = 1:T
    configstar(i) = treeMAX{i}.alphabet{treeMAX{i}.xstar};
end
configstar
MAP = treeMAX{1}.msgs{1}



% run sum-product algorithm
treeSUM = sumprod(tree);

% print the results
strHID = '%4.2f %6.2f %6.2f\n';
strVIS = '%4.2f %6.2f\n';
for i = 1:length(treeSUM)
    [y ind] = max(exp(treeSUM{i}.logprobs));
    if i>T; str = strVIS; else str = strHID; end
    fprintf(str,exp(treeSUM{i}.logprobs));
    foo(i) = treeSUM{i}.alphabet{ind};
end


foo
observations
