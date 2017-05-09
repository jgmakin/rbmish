clear; clc; close all
path(path,'..\glmnet');

monkey = 'E';

% init
day = 1;
options = glmnetSet;
alpha = 0.01:0.01:0.2;

% loop through days and data
while 1 % for k = 1:1 % 
    for trgnum = 1:3
        
        [X y angle] = dataExtractor2(day,trgnum,monkey);
        if angle == inf
            break; 
        else
            RR = [];
            angleCell{day}{trgnum} = angle;
            % loop through alphas
            for i = 1:length(alpha)
                % a = 0.03;
                options.alpha = alpha(i);
                
                % fit = glmnet(X,y,'gaussian',options);
                % results
                % glmnetPrint(fit);
                % glmnetPlot(fit,'dev','true');
                % glmnetCoef(fit,0.01);
                % glmnetPredict(fit,'response',X(1:10,:),[0.01,0.005]') % make predictions
                
                % figure()
                
                
                % elastic net w/cross-validatation
                fprintf('angle: %d; ',angle);
                CVerr = cvglmnet(X,y,size(X,1),[],'response','gaussian',options,0);
                RR{i} = CVerr.Rsquared';
            end
            
            RRR{trgnum} = RR;
            fprintf('\n\n');
        end
    end
    
    if angle == inf 
        break;
    else
        RRRR{day} = RRR;
        day = day + 1;
    end
    
end

% save('cvRsqDmitri.mat','RRRR','angleCell','alpha ');




% RRRR: {day}{trgnum}{alpha}{#of nz beta components}
k = 1;
for day = 1:length(RRRR)
    for trgnum = 1:length(RRRR{day})
        for i_alpha = 1:length(RRRR{day}{trgnum})
            v(i_alpha) = max(RRRR{day}{trgnum}{i_alpha});
        end


       
        if 0
            bar(alpha,v);
            titlestr = ['angle: ',num2str(angleCell{day}{trgnum})];
            titlestr = [titlestr ', day: ',num2str(day)];
            titlestr = [titlestr ', trgnum: ',num2str(trgnum)];
            axis([0, alpha(end)+alpha(1), min(v), max(v)]);
            title(titlestr);
            pause()
        end
        
      %  bar(alpha,v);
      %  axis([0, alpha(end)+alpha(1), min(v), max(v)]);
     %  pause()    
            
%         if trgnum == 2
%            cntrtrg{k} = v;
%         end
        peakRsq(k) = max(v);
        angleVec(k) = angleCell{day}{trgnum};
        k = k+1;
    end
    
end



% close
figure
scatter(angleVec,peakRsq);
hold on
ANG = reshape(angleVec,3,length(angleVec)/3);
ANG(3,:) = ANG(3,:) + (ANG(1,:) == 270)*360;        % 0 -> 360 for some 0s
PRS = reshape(peakRsq,3,length(peakRsq)/3);
plot(ANG,PRS);
hold off;

% 315       7:1         0.03
% 225       7:2         0.17















