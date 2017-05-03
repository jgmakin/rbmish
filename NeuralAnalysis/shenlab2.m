% sparse coding of 


clc; close all;
path(path,'C:\#code\sparsenet');


% extract data
[U,X,angle] = dataExtractor2(5,1,'E');

% set params
[Nsamples Nneurons] = size(U);
Nbehavior = size(X,2);
Ncoeff = 20;

% sparse coding (assume x = U*Q*s, w/x behavioral data and U neural data)
Q = rand(Nneurons, Ncoeff) - 0.5;
A = U*Q;
A = A*diag(1./sqrt(sum(A.*A)));
figure(1); colormap(gray);
figure(2); colormap(gray);
sparsenetMod