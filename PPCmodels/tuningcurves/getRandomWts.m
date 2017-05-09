function wts = getRandomWts(wts,Nhid,Nvis)
% randomly populate the weights, but with the same mean and variance as the
% learned feedforward weights (altho' these are leptokurtotic).  The
% "random" biases are just initialized at zeros.

%-------------------------------------------------------------------------%
% Created: 07/27/12
%   by JGM
%-------------------------------------------------------------------------%


% from the loaded wts
vishid = wts{1}(1:end-1,:);

% statistics of the feedforward weights
std = sqrt(var(vishid(:)));
mu = mean(vishid(:));

% random wts
vishid      = mu + std*randn(Nvis, Nhid);
hidbiases   = 0*ones(1,Nhid);
visbiases   = 0*ones(Nvis,1);
wts{1} = [vishid; hidbiases];
wts{2} = [vishid'; visbiases'];

end