function plotTCs(s,V,params)
% plot 2-dimensional tuning curves on MxM grid, with (hidden-unit) data V.

%-------------------------------------------------------------------------%
% Created: 07/27/12
%   by JGM
%-------------------------------------------------------------------------%


% get dimensionality of the stimulus
Ndims = params.Ndims;

% params
n = 8;
clrs = 'rgbmcy';

% sort the x-axis data (stimuli) for plotting
% [foo sorta] = sort(x);

% get some random indices for plotting
[y,ind] = sort(rand(1,size(V,2)));
% ind = 1:size(V,2);
%%%%%%%%
ind = 1:n^2;
%%%%%%%%

figure(4);
for j = 1:size(V,3)
    
    for i = 1:n^2
        
        subplot(n,n,i);
        if Ndims==1
            hold on;
            plot(s,V(:,ind(i),j),clrs(j)); 
            axis([-5 5 0 1]);
            hold off;
        else
            % put back into long format---in this case, gridded by theta
            M = sqrt(size(V,1));
            vv = shortdata(M,3,V(:,:,j));
            xx = [s(1,1),s(end,1)];
            yy = [s(1,2),s(end,2)];
            imagesc(xx,yy,squeeze(vv(:,ind(i),:)),[0,1]);
            axis xy; axis tight;
           % axis([params.thmin(1) params.thmax(1) params.thmin(2) params.thmax(2)]);
        end       
    end
    
    if (Ndims==2)&&(j<size(V,3)), figure; end % pause(); end;
end


end