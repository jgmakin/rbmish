function EdgeMat = getEdgeMatrix(start,finish,dt,m)
% get a matrix of edges, corresponding to a sliding window, for histc.

%-------------------------------------------------------------------------%
% Created: 05/14/13
%   by JGM
%-------------------------------------------------------------------------%

binsize = m*dt;
Nwindows = length(start:binsize:finish);

edge0 = start - binsize/2;
edgeF = edge0+Nwindows*binsize; % min((edge0+Nwindows*binsize),finish);

EdgeMat = repmat([edge0:binsize:edgeF],m,1) + (1:m)'*ones(1,Nwindows+1)*dt;
EdgeMat(EdgeMat>finish) = finish;
EdgeMat(EdgeMat<start) = start;


if 0
    foo = colormap;
    for offset=1:m
        figure(166); hold on;
        scatter(EdgeMat(offset,:),offset*ones(1,Nwindows+1),[],foo(offset*4,:));
    end
    pause()
end

end




% foo = colormap;
% figure(166); clf; hold on;
% 
% for offset=1:m
%     
%     
%     figure(166); hold on;
%     scatter(EdgeMat(offset,:),offset*ones(1,Nwindows+1),[],foo(offset*4,:));
%     
%     
% end
% hold off;