function [X,Ns] = plotEMMfit(Y,EMMparams,Nbins,fignum)
% Plot Erlang mixture model data (histogram) and distribution (line)
% 
% USAGE:
%   [Pizq,Thyq] = getEMMposteriorProbs(testData.R,wts,params);
%   EMMparams = EFparams2ErlangParams(Pizq,Thyq);
%   [X,Ns] = plotEMMfit(Y,EMMparams,100,fignum);

%-------------------------------------------------------------------------%
% Cribbed: 06/03/16
%   from testEFHdecoding
%   by JGM
%-------------------------------------------------------------------------%

% histogram of data
figure(fignum); clf; hold on;
X = linspace(0,1,Nbins);
Ns = hist(Y(:),X);
Ns = (Ns./sum(Ns,2))'/(X(2)-X(1));    % normalize
bar(X,Ns,'FaceColor',[253,180,98]/255,'EdgeColor','none');

% probability density function
if ~isempty(EMMparams)
    x = 0:0.005:1;
    M = length(x);
    y = arrayfun(@(ii)(EMMparams(ii).pis'*gampdf(...
        repmat(x,[length(EMMparams(ii).pis),1]),...
        repmat(EMMparams(ii).ks,[1,M]),...
        repmat(EMMparams(ii).mus,[1,M]))),...
        1:length(EMMparams),'UniformOutput',false);
    plot(x,cat(1,y{:}));
end

% figure props
xlabel('coherence');
ylabel('probability density');
axis tight;
hold off;

end

