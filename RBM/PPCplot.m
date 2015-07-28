function PPCplot(plotdata,params,titlestring)
% PPCPLOT    Population plots
%   PPCPLOT plots PPCs in an appropriate way

%-------------------------------------------------------------------------%
% Created: 1/13/10
%   by JGM
%-------------------------------------------------------------------------%

switch params.Ndims
    case 1
        scatter(1:length(plotdata(:)),plotdata(:)); % axis off;
    case 2
        imagesc(plotdata); axis equal; axis off;
    otherwise
        fprintf('you never worked out PPCplots for m!=2 -- jgm\n');
end

title(titlestring);

end


%%% for making Gaussians --- insert into respfxn.m
% x = repmat(b,params.N,1);
% y = repmat(b',params.N,1);
% y = y(:);
% % y = y(length(y):-1:1);
% 
% 
% 
% 
% scatter3(y,x,r,'k'); 
% hold on;
% colormap gray;
% % surfl(b,b,g*f); shading interp; colormap gray;
% surf(b,b,g*f);
% grid off; axis off;
% hold off;
% 5