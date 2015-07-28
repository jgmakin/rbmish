function [aopt bopt] = templatematcher(noisyhill,noiselesshill,METRIC)
% TEMPLATEMATCHER   Fits bubbles to bubbles
%   The maximum-likelihood estimate for a population of independent units
%   is given by "the noise-free hill that minimizes distance from the data,
%   where the distance metric is determined by the distribution of the
%   noise...."
%   "[I]f the noise is Gaussian, the appropriate distance is the
%   Mahalanobis norm...."  --And what about Poisson noise??


N = size(noisyhill,1);

aopt = 0; bopt = 0;
d = distance(zeros(size(noisyhill)),noisyhill,METRIC);

%%% NOTE %%%
% You can probably orthogonalize this; i.e., do each direction
% independently.....right?
%%%%%%%%%%%%
for i = 1:N
    for j = 1:N
        % input
        fmat = noiselesshill(1+i-1:N+i-1,1+j-1:N+j-1);
        dnew = distance(fmat,noisyhill,METRIC);
        if dnew < d
            d = dnew;
            aopt = i; bopt = j;
        end
    end
end


end


function d = distance(A,B,METRIC)

if strcmp(METRIC,'Euclidean');
    d = sqrt(sum(sum((A - B).^2)));
    % elseif
else
    error('unrecognized distance metric -- jgm');
end

end