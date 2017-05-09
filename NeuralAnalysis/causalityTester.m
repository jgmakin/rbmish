function causalityTester(TD)

clc
fprintf('"tuning series"\n');

for i = 1:length(TD)
    
    Ntrials = size(TD(i).tunRATES,1);
    y = cosd(TD(i).tunIA);
    X = [ones(Ntrials,1) TD(i).tunRATES]; 
    [beta yhat RsqYX] = regressPCA(y,X);
    
    res = y - yhat;
    % W = [ones(Ntrials,1) res];
    W = [ones(Ntrials,1) res res.^2 res.^3 res.^4];
    z = cosd(TD(i).tunTRG);
    [alpha zhat RsqZW] = regressPCA(z,W);
    
    fprintf('R^2''s: IA on rates: %f; TRG on residuals: %f\n',RsqYX,RsqZW);
    
%     alpha = w'*y\w'*z;
%     zhat = w*alpha;
%     res = z - zhat;
%     
%     SSerr = res'*res;
%     zcntr = z - mean(zhat);
%     SStot = zcntr'*zcntr;
%     Rsq = 1 - SSerr/SStot;

    
end

fprintf('\n\n"experimental series"\n');


for i = 1:length(TD)
    
    Ntrials = size(TD(i).expRATES,1);
    y = cosd(TD(i).expIA);
    X = [ones(Ntrials,1) TD(i).expRATES]; 
    [beta yhat RsqYX] = regressPCA(y,X);
    
    res = y - yhat;
    % W = [ones(Ntrials,1) res];
    W = [ones(Ntrials,1) res res.^2 res.^3 res.^4];
    z = cosd(TD(i).expTRG);
    [alpha zhat RsqZW] = regressPCA(z,W);
    
    fprintf('R^2''s: IA on rates: %f; TRG on residuals: %f\n',RsqYX,RsqZW);
    
%     alpha = w'*y\w'*z;
%     zhat = w*alpha;
%     res = z - zhat;
%     
%     SSerr = res'*res;
%     zcntr = z - mean(zhat);
%     SStot = zcntr'*zcntr;
%     Rsq = 1 - SSerr/SStot;

    
end



end


% (1) cross-validation
% (2) fix experimental trials!