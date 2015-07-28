function EFHdisp(DISP,pvisstates,qhidstates,qvisstates,...
    vishid,rows,cols,indices,params,ax)
% EFHdisp   EFH display
%   EFHdisp displays useful information associated with the currently
%   training EFH.  This file was created just to unclutter EFH.m

%-------------------------------------------------------------------------%
% Revised: 06/03/15
%   -modernized, somewhat
% Revised: 09/18/14
%   -renamed rbmdisp.m -> EFHdisp.m 
% Cribbed: 01/28/11
%   -from rbm.m
%   by JGM
%-------------------------------------------------------------------------%

if DISP(1)
    % display visible data and (visible) confabulations (for case 1)
    figure(2); colormap(gray);
    Ttrue = displayshape(pvisstates(1,:),params);
    subplot(2,1,1);
    PPCplot(cat(2,Ttrue{:}),params,'input responses');
    Tconfab = displayshape(qvisstates(1,:),params);
    subplot(2,1,2);
    PPCplot(cat(2,Tconfab{:}),params,'output responses');
end

if DISP(2)
    % display hidden states
    figure(3); colormap(gray);
    imagesc(qhidstates)
end

if DISP(3)
    % display weights
    for ii = 1:rows
        for jj = 1:cols
            Twts = displayshape(vishid(:,indices((ii-1)*rows+jj)),params);
            imagesc([Twts{1},Twts{2}],'Parent',ax(ii,jj));
            set(ax(ii,jj),'xtick',[],'ytick',[]);
        end
        
    end
end
% title('weights');

% if DISP(4)
%     fprintf('max <W_i>: %f, min <W_i>: %f\n',...
%         max(mean(vishid,2)),min(mean(vishid,2)));
%     fprintf('max{hb}: %f, min{hb}: %f\n',...
%         max(hidbiases),min(hidbiases));
%     fprintf('max{vb}: %f, min{vb}: %f\n',...
%         max(visbiases),min(visbiases));
%     fprintf('max{d_c}: %f, min{d_c}: %f\n',...
%         max(max((negvisSfctStats))),min(min((qvisstates))));
%     fprintf('max{sum(d)}: %f, max{sum(d_c)}: %f\n',...
%         max(posvisact),max(negvisact));
%     fprintf('max{h*d}: %f\n',max(max(negprods)));
% end

drawnow;

end