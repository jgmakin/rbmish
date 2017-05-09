function net = nnFit(inputs,targets)
%CREATE_FIT_NET Creates and trains a fitting neural network.
%
%  NET = CREATE_FIT_NET(INPUTS,TARGETS) takes these arguments:
%    INPUTS - RxQ matrix of Q R-element input samples
%    TARGETS - SxQ matrix of Q S-element associated target samples
%  arranged as columns, and returns these results:
%    NET - The trained neural network
%
%  For example, to solve the Simple Fit dataset problem with this function:
%
%    load simplefit_dataset
%    net = create_fit_net(simplefitInputs,simplefitTargets);
%    simplefitOutputs = sim(net,simplefitInputs);
%
%  To reproduce the results you obtained in NFTOOL:
%
%    net = create_fit_net(Y,th);

% Create Network
% numHiddenNeurons = 20;  % Adjust as desired
% net = fitnet([10 5]);
net = fitnet([20 15]);
net.divideParam.trainRatio = 70/100;  % Adjust as desired
net.divideParam.valRatio = 15/100;  % Adjust as desired
net.divideParam.testRatio = 15/100;  % Adjust as desired

% Train and Apply Network
[net,tr] = train(net,inputs,targets);
outputs = net(inputs);

% Plot
plotperf(tr)
plotfit(net,inputs,targets)
plotregression(targets,outputs)
