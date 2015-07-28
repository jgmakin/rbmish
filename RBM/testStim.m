function testStim(D,X,wts,params)
% UNIMODALSTIM  Test unimodal input to a trained DBN
%  USAGE: 
%       testStim(testbatchdata0,xtest,wts,params);
%       testStim(batchdata0,xtrain,wts,params);
%  UNIMODALSTIM takes a dataset D (either batchdata0 or testbatchdata0) and
%  the *corresponding* causal stimuli (xtrain or xtest, resp.) along with
%  the weights of a trained-up DBN and its params, and generates images and
%  statistics of the inputs and outputs.
%-------------------------------------------------------------------------%
% Revised: 08/03/10
%   -add params as an argument to updown
% Created: 07/13/10
%   by JGM
%-------------------------------------------------------------------------%

% params
ALTINPUT = 'Hand-Position';                               % vis/prop;
ALT = 'degrade';                                  % shift/kill/degrade;

% wipe out one input
i = strcmp(ALTINPUT,'Joint-Angle');
if strcmp(ALT,'kill')
    data = killinput(D,i);
elseif strcmp(ALT,'shift')
    data = shiftinput(D,i,params);
elseif strcmp(ALT,'degrade')
    data = degradeinput(D,i);
else
    error('unrecognized input alteration--- jgm');
end
    
% push it through the model (in mean mode) 
[tilde,dataout] = updownfast(data,wts,params);

% display the results
griddisp(data(:,:,end),dataout,1000,params);

% compute some statistics
[covVi,covPi,covVo covPo,mu] = estimatorStats(X,data,wts,params)

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function halfdata = killinput(D,m)

halfdata = zeros(size(D));
for i=1:size(D,1)
    for j = 1:size(D,3)
        T = displayshape(D(i,:,j),params);
        % fooP = zeros(size(fooP));
        halfdata(i,:,j) = [T{1}(:)*m; T{2}(:)*(1-m)]';
    end
end


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function shiftdata = shiftinput(D,m,params)

% init
r = 2*params.c;                                 % radius as fxn(std. dev.)
angle = rand*2*pi;                              % and in some random direc.
gran = params.granularity;
shift = round(gran*r*[cos(angle) sin(angle)]);

% shift
shiftdata = zeros(size(D));
for i=1:size(D,1)
    for j = 1:size(D,3)
        T = displayshape(D(i,:,j),params);
        T{1} = circshift(T{1},(1-m)*shift);
        T{2} = circshift(T{2},m*shift);
        shiftdata(i,:,j) = [T{1}(:); T{2}(:)]';
    end
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function degradeddata = degradeinput(D,m)

degradeddata = zeros(size(D));
e = 0.1*randn(size(D,2)/2,1);                   % recall: data \in [0,1]

for i=1:size(D,1)
    for j = 1:size(D,3)
        T = displayshape(D(i,:,j),params);
        degradeddata(i,:,j) = [T{1}(:) + e*(1-m); T{2}(:) + e*m]';
    end
end


end
%-------------------------------------------------------------------------%