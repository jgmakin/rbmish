function [R,Q] = getDataBouncingBalls(X,Q,balls,res)
% makeBouncingBalls     Generate time series of 2D ricocheting balls
%
% USAGE:
%   [R,Q] = getDataBouncingBalls(X,Q,balls,resolution);

%-------------------------------------------------------------------------%
% Adapted: 08/10/15
%   by JGM
%   -from Ilya Sutskever's python code
%-------------------------------------------------------------------------%

% only "encode" positions
X = reshape(X,size(X,1),[],balls.N);
X = X(:,1:2,:); 

% Ns
[Nexamples,Ndims,Nballs] = size(X);

% default params
if ~isfield(balls,'r'), balls.r = ones(1,Nballs)*1.2; end
if ~isfield(balls,'boundingbox'), balls.boundingbox = 10; end

% (Nexamples x Ndims x Nballs) -> (1 x 1 x Nexamples x Nballs x Ndims)
X = permute(X,[5,4,1,3,2]);

% embed into Gaussian function
[I,J] = meshgrid(ar(0,1,1./res)*balls.boundingbox,...
    ar(0,1,1./res)*balls.boundingbox);
r = reshape(balls.r,[1,1,1,Nballs]);
A = sum(exp(-(((I - X(:,:,:,:,1)).^2 + (J - X(:,:,:,:,2)).^2).*r.^2).^4),4);

% saturate
A(A>1)=1;

R = reshape(A,[res^Ndims,Nexamples])';

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function foo = ar(x,y,z)

foo = z/2 + x:z:y;

end
%-------------------------------------------------------------------------%


% %-------------------------------------------------------------------------%
% function A = bounce_mat(res, n, T, r, m)
% if isempty(r), r = ones(n,1)*1.2; end
% if isempty(n), n = 2; end
% if isempty(T), T = 128; end
% 
% SIZE = 10;
% x = bounce_n(T,n,r,m,SIZE);
% A = matricize(x,res,r);
% end
% %-------------------------------------------------------------------------%
% 
% 
% %-------------------------------------------------------------------------%
% function vec = bounce_vec(res, n, T, r, m)
% if isempty(r), r = ones(n,1)*1.2; end
% if isempty(n), n = 2; end
% if isempty(T), T = 128; end
% 
% SIZE = 10;
% x = bounce_n(T,n,r,m,SIZE);
% V = matricize(x,res,r);
% vec = reshape(V,[T, res^2]);
% end
% %-------------------------------------------------------------------------%


% %-------------------------------------------------------------------------%
% function show_single_V(V)
% res = round(sqrt(shape(V)[0])) %%%%
% show(V.reshape(res, res)); %%%%%
% end
% %-------------------------------------------------------------------------%

