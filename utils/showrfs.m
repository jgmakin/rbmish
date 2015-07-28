function showrfs(A,bg)
% Plot 2D receptive fields in the usual way.
%   showrfs expects a matrix A whose *columns* are examples---i.e.,
%   neurons.  Each column should have a perfect-square number of elements
%   (i.e., the number of rows of A), since the receptive field was
%   presumably mapped out at an equal number of discrete points in x and y.

%-------------------------------------------------------------------------%
% Revised: 11/04/13
%   -added axis xy to flip the image into "normal" coordinates!!
% Cribbed: 11/04/13
%   -from code by B.A. Olshausen (c. 2006)
%-------------------------------------------------------------------------%



if ~exist('bg','var')
    bg='black';
end

[L,M]=size(A);

sz=sqrt(L);

buf=1;

if floor(sqrt(M))^2 ~= M
    %%% m=sqrt(M/2);
    m = floor(sqrt(M));
    n=M/m;
    %%% of course this won't always work!
else
    m=sqrt(M);
    n=m;
end

if strcmp(bg,'black')
    array=-ones(buf+m*(sz+buf),buf+n*(sz+buf));
else
    array=ones(buf+m*(sz+buf),buf+n*(sz+buf));
end


k=1;
for j=1:n
    for i=1:m
        
        clim=max(abs(A(:,k)));
        array(buf+(i-1)*(sz+buf)+[1:sz],buf+(j-1)*(sz+buf)+[1:sz])=...
            reshape(A(:,k),sz,sz)/clim;
        k=k+1;
        
    end
end

imagesc(array,[0 1])
axis image off
axis xy
colormap('gray')

end
%-------------------------------------------------------------------------%