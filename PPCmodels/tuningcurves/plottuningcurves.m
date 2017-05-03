function plottuningcurves(wts,params)
% explore the hidden unit tuning curves.  What do they look like?  Do they
% show "inverse effectiveness"? sub- or super-additivity?  Etc.
%
% NB: This was more or less hard-coded for params.Ndims=1, length(params.mods)=3.
%
% First run:
%   load results/new/add021412.mat
% then call
%   plottuningcurves(wts,params)

%-------------------------------------------------------------------------%
% Revised: 12/16/13
%   -X -> S and corresponding changes
% Revised: 07/20/12
%   -replaced old delta-finding fxn w/new one ("getDeltaDstrb2"), etc.
% Revised: 07/18/12
%   -functionized, etc.
% Created: 07/17/12
%   by JGM
%-------------------------------------------------------------------------%


% params
Nexamples = 1000;
hidDstrbs = params.typeUnits{2};
hidNums = params.numsUnits{2};
Nhid = sum(hidNums);
USERANDOMWTS = 0;
gazes = (1:3)/4;


% use random wts?
if USERANDOMWTS
    Nvis = sum(params.numsUnits{1});
    wts = getRandomWts(wts,Nhid,Nvis);
end

% generate a standard set of stimuli to play with
[~,S0] = generateData(Nexamples,params); % ,'swing',0);

% sort the data by ascending hand position (x)
[~,indmax] = sort(S(:,:,strcmp(params.mods,'Hand-Position')));
S = S(indmax,:,:);


% generate some noiseless data with a fixed eye position

% malloc 
V = zeros(size(S,1),Nhid,length(gazes));
gains = mean([params.gmin; params.gmax]);
params.gmin = gains;
params.gmax = gains;
%%% shall I fix the gains??
for iGaze = 1:length(gazes)
    
    S = adjustForNewGaze(S,gazes(iGaze),params);
    params.typeUnits{1} = {'Dirac'};
    R = generateData(Nexamples,params,'stimuli',S);
    V(:,:,iGaze) = invParamMap(R,wts{1}(1:end-1,:),wts{1}(end,:),...
        hidDstrbs,hidNums,params);
end

% plot
plotTCs(S(:,:,1),V,params);

% get distribution of LMMM's "deltas"
% getDeltaDstrb(gazes,[1,3],X(:,1),V,params);

% try a "template matching" kind of approach...
[deltas,shiftable] = getDeltaDstrb2(gazes,[1,3],S(:,:,1),V,params);

% plot the incoming weights to each hidden unit
plotInWts(deltas,shiftable,wts,params)

end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%




%-------------------------------------------------------------------------%
function S = adjustForNewGaze(S,eyepos,params)
% keeping vis (S(:,:,1)) the same, change the old gaze angle(s) (S(:,:,3)) 
% to the new one, eyepos, and then adjust prop (S(:,:,2)) accordingly.

% init
[Nexamples, Ndims, Nmods] = size(S);
if Nmods~=3, error('wrong number of modalities -- jgm'); end
eyemin = params.roboparams.eyemin; 
eyemax = params.roboparams.eyemax;

% replace old eyeposition with new one
e = scalefxn(eyepos,zeros(size(eyemin)),ones(size(eyemax)),eyemin,eyemax);
S(:,:,strcmp(params.mods,'Gaze-Angle')) = repmat(e,[Nexamples,1,1]);

% correct theta accordingly
S(:,:,strcmp(params.mods,'Joint-Angle')) =...
    IK2link(S(:,:,strcmp(params.mods,'Hand-Position')) -...
    S(:,:,strcmp(params.mods,'Gaze-Angle')),params.roboparams,1);
%%% changed but never checked: 7/2/14 jgm

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function deltas = getDeltaDstrb(allgazes,gazeinds,x,V,params)

% params
thresh = 0.15;              % throw out tuning curves that never reach this
Ngazes = 2;
Nhid = size(V,2);
eyemin = params.roboparams.eyemin; 
eyemax = params.roboparams.eyemax;
Ndims = params.Ndims;


% malloc
argmax = zeros(Ngazes,Nhid);
negligible = zeros(Nhid,1);
eye = zeros(Ngazes,1);

% loop through hidden units and gazes
for hiddenunit = 1:Nhid
    for gaze = 1:Ngazes
        [fxnmax,indmax] = max(V(:,hiddenunit,gazeinds(gaze)));
        argmax(gaze,hiddenunit) = x(indmax);
        if (fxnmax < thresh), negligible(hiddenunit) = 1; end;
        
        eye(gaze) = scalefxn(allgazes(gazeinds(gaze))*ones(1,length(eyemin)),...
            zeros(size(eyemin)),ones(size(eyemin)),eyemin,eyemax);
    end
end

% compute "deltas"
eyediff = diff(eye);

deltas = diff(argmax,1,2);                      % throw out low TCs and
deltas = deltas(~negligible)/eyediff;           %   scale into 0-1 range

% histogram
hist(deltas,20);


end
%------------------------------------------------------------------------%


%------------------------------------------------------------------------%
function [deltas,shiftable] = getDeltaDstrb2(allgazes,gazeinds,x,V,params)


% params
[Nexamples,Nhid,Nshifts] = size(V);
if length(allgazes)~=Nshifts, error('gazeinds vector wrong size -- jgm'); end
thresh = 0.15;              % throw out tuning curves that never reach this


% malloc
bestIndShfts = zeros(Nhid,1);
bestMSEs = zeros(Nhid,1);

shifts = -(Nexamples-1):(Nexamples-1);
pad = 50;
tic
for hiddenunit = 1:Nhid
    
    % tuning curves
    v1 = V(:,hiddenunit,gazeinds(1));
    v2 = V(:,hiddenunit,gazeinds(2));

    if (max(v1) > thresh) && (max(v2) > thresh)
        
        bestMSE = Inf;
        % bestCrln = 0;
        for i = (1+pad):(length(shifts)-pad)% 1:length(shifts)
            
            % some clever indexing
            z = shifts(i); 
            inds2 = (1+max(z,0)):(Nexamples+min(z,0));
            inds1 = fliplr(Nexamples+1-inds2);
            
            % compute MSE for this shift
            xx = v1(inds1); yy = v2(inds2);
            MSE = mean((xx-yy).^2);
            
            % crln = (xx'*yy)/(sqrt(xx'*xx)*sqrt(yy'*yy));
            
            % if it's better than its best predecessor, it's the best
            if (MSE<bestMSE) && (max(v1(inds1))>0.15) && (max(v2(inds2))>0.15)
            % if (crln>bestCrln) && (max(v1(inds1))>0.15) && (max(v2(inds2))>0.15)
                bestMSE = MSE;
                bestShft = z;
            end
        end
        bestIndShfts(hiddenunit) = bestShft;
        bestMSEs(hiddenunit) = bestMSE;
    end
        
end     
toc
        

% plot
shfts = plotShiftedTuningCurves(bestIndShfts,shifts,x,V,gazeinds);


% histogram the deltas
deltas = deltahist(shfts,allgazes(gazeinds),params);

% indices
shiftable = find(bestIndShfts);

end
%------------------------------------------------------------------------%


%------------------------------------------------------------------------%
function shfts = plotShiftedTuningCurves(indShfts,shifts,x,V,gazeinds)

% find the shifts for tuned neurons
inds = find(indShfts~=0);

% find the corresponding tuning curves
V1 = V(:,inds,gazeinds(1)); %%% should correspond to gaze angle at 1/4
V3 = V(:,inds,gazeinds(2)); %%% should correspond to gaze angle at 3/4

% transform the index shifts into cm (or whatever) shifts
len = max(x) - min(x);
oo = ones(length(inds),1);
shfts = scalefxn(indShfts(inds),min(shifts)*oo,max(shifts)*oo,-len*oo,+len*oo);

% plot the shifted and unshifted tuning curves
figure(48)
for i = 1:min(length(inds),49) 
    subplot(7,7,i), 
    hold on, 
    plot(x,V1(:,i),'b'); 
    plot(x,V3(:,i),'r');
    plot(x+shfts(i),V1(:,i),'g');
end



end
%------------------------------------------------------------------------%


%------------------------------------------------------------------------%
function deltas = deltahist(shfts,gazes,params)

% params
numbins = 20;
eyemin = params.roboparams.eyemin;
eyemax = params.roboparams.eyemax;
oo = ones(size(eyemin));

% gaze angles in cm (or whatever)
eye1 = scalefxn(gazes(1)*oo(:)',0*oo,1*00,eyemin,eyemax);
eye2 = scalefxn(gazes(2)*oo(:)',0*oo,1*00,eyemin,eyemax);
eyediff = eye2 - eye1;
deltas = shfts/eyediff;

figure;
hist(deltas,numbins);

figure(48);
for i = 1:min(length(deltas),49)
    subplot(7,7,i);
    title(deltas(i));
end

end
%------------------------------------------------------------------------%


%------------------------------------------------------------------------%
function plotInWts(deltas,shiftable,wts,params)

% init
Nmods = length(params.mods);
N = params.N;
Nhid = size(wts{1},2);
thresh = 0.25;

% make some figures
h = zeros(4,1);
for i = 1:4, h(i) = figure(i+76);end

% cycle through all the hidden units
modwt = zeros(length(params.mods),1);
n = 10;
for i = 1:Nhid
    
    % what color?
    if any(shiftable==i)
        if abs(deltas(shiftable==i)-1) < thresh
            clr = 'g'; 
        elseif abs(deltas(shiftable==i)) < thresh
            clr = 'm';
        else 
            clr = 'r';        
        end
    else
        clr = 'b';
    end
    
    
    % the incoming weight vector to this hidden unit
    w = wts{1}(:,i);

    if i>n^2, f = h(2); else f = h(1); end
    set(0,'CurrentFigure',f);
    subplot(n,n,mod(i-1,100)+1);      
    plot(w,clr);
    
    for j = 1:Nmods
        vec = (1:N) + (j-1)*N;
        modwt(j) = sum(abs(w(vec))); 
    end
         
    if i>n^2, f = h(4); else f = h(3); end
    set(0,'CurrentFigure',f);
    subplot(n,n,mod(i-1,100)+1);      
    bar(modwt,clr);
    
end

end
%------------------------------------------------------------------------%



%     % get the relevant tuning curves and zero pad
%     v1 = [V(:,hiddenunit,1); zeros(length(shifts),1)];
%     v2 = [V(:,hiddenunit,3); zeros(length(shifts),1)];
%         
           
        
%         % try all the shifts of v1
%         row = [fliplr(v1(1:nExamples)'), zeros(1,nExamples-1)];
%         col = circshift(v1,shifts(1));
%         C1 = toeplitz(col,row);
%         C2 = repmat(v2,1,size(C1,2));
%         M = C2 - C1;
%         M = M(1:(2*nExamples-1),:);
%         [MSE1 ind1] = min(sum(M.^2,1));
%         
%         % now try all the shifts of v2
%         row = [fliplr(v2(1:nExamples)'), zeros(1,nExamples-1)];
%         col = circshift(v2,shifts(1));
%         C2 = toeplitz(col,row);
%         C1 = repmat(v1,1,size(C1,2));
%         M = C2 - C1;
%         M = M(1:(2*nExamples-1),:);
%         [MSE2 ind2] = min(sum(M.^2,1));
%         %%%% need to scale this by the number of data...altho' it works
%         %%%% ok w/o.....
        
        % use which ever was better
%         delta(hiddenunit) = (MSE1<=MSE2)*shifts(ind1)+(MSE1>MSE2)*shifts(ind2);
%         LorR(hiddenunit) = (MSE1<=MSE2) - (MSE1>MSE2);
    % end 

