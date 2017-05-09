function isomap4PNSdata(S,UnitSpikes,fignum)
% isomap4PNSdata    Run the Kiani2015 analysis on data from the Sabes lab
%
% USAGE:
%{
    datadir = 'C:\#code\HHS\extracteddata\';
    tag = 'D080620';
    load([datadir,'KFtuningdataHHS',tag]); % "the tuning series"
    isomap4PNSdata(St,UnitSpikesT,8888)
%}

%-------------------------------------------------------------------------%
% Created: 03/27/15
%   by JGM
%-------------------------------------------------------------------------%


% useful params
KFparams.Nstates = 6;
KFparams.Ndims = 2;
KFparams.m = 16;                              % 16 => 60 Hz (66.7 ms bins)
KFparams.dt = 1/240;                          %
KFparams.BINMETHOD = 'nonoverlappingwindow';
numNbs = 1;
kmeansk = 2;


% bin the spike counts and resample the trajectories
[R,Xout,endinds] = binSpikeCounts(S,UnitSpikes,KFparams);
vec = mean(R)>0;
R = R(:,vec);
unitIDs = cat(1,UnitSpikes(:).id);
unitIDs = unitIDs(vec,:);

% isomap, k-means
d = getDissimilarityMatrix(R);
options.dims = 1:10;


% keep increasing the number of neighbors till isomap finds one clique
NmaxClique = 0;
while NmaxClique~=sum(vec)
    
    % start by trying 2
    numNbs = numNbs+1;
    [Y,Res,E] = Isomap(d,'k',numNbs,options,fignum);
    NmaxClique = length(Y.index);
end



% now put back onto the cortex---using the k groups from k means...
[idx,ctrs] = kmeansFitAndPlot(Y.coords{3}',kmeansk,fignum);
colorElectrodeArrayByGroup(idx,unitIDs,fignum);

% ...or using a polar color map...
C = smoothlyColor(Y.coords{3}',fignum);
colorElectrodeArrayByPolarCoords(C,unitIDs,fignum)






end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function d = getDissimilarityMatrix(Y)
% turn array firing rates into a "dissimilarity matrix"

% (pairwise) correlation of HHS neural data
Ycntrd = Y - mean(Y,1);
Ystd = sqrt(sum(Ycntrd.^2));
Yzscored = Ycntrd./Ystd;

% weird, naive-seeming "distance" measure of Kiani et alia
Ycorr = Yzscored'*Yzscored;
d = 1 - Ycorr;

% alternatively, you might use the matrix of distances
% dSq = bsxfun(sum(Yzscored.*Yzscored) + sum(Yzscored.*Yzscored)') - 2*(Yzscored'*Yzscored);
% dSq(dSq<0) = 0;
% d = sqrt(dSq);
%%% the results are *very* similar; but actually this is slightly worse


end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function C = smoothlyColor(X,fignum)


% asign a color to each point based on polar coordinates
rmax = max(sqrt(sum(X(:,1:2).^2,2)));
C = [(atan2(X(:,2),X(:,1))+pi)/2/pi,....
    sqrt(sum(X(:,1:2).^2,2))/rmax,...
    ones(size(X,1),1)];
C = hsv2rgb(C);


% plot
figure(fignum);
subplot(2,3,6);
scatter(X(:,1),X(:,2),[],C(:,:));



end
%-------------------------------------------------------------------------%



%-------------------------------------------------------------------------%
function [idx,ctrs] = kmeansFitAndPlot(X,k,fignum)
% just like it says

% colors!  You never expect to need more than six
clrs = 'rgbymk';

% fit with k-means
[idx,ctrs] = kmeans(X, k);

% plot---in 2D in any case
figure(fignum); 
subplot(2,3,5); cla;
hold on;
for iGroup = 1:k
    scatter(X(idx==iGroup,1), X(idx==iGroup,2),clrs(iGroup));
end
scatter(ctrs(:,1),ctrs(:,2),'kx');
hold off;

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function colorElectrodeArrayByPolarCoords(C,unitIDs,fignum)
% color the electrode array according to the cluster ids


% from the spreadsheet data for Dmitri's 3/24/08 PMd array (1024-0217)
pin2elec = [93,92,94,95,75,96,85,97,86,98,87,88,77,99,66,89,76,90,67,79,...
    58,80,78,70,68,60,69,50,59,40,49,30,83,84,73,74,63,64,53,54,43,55,...
    44,45,33,46,34,65,24,56,35,47,25,57,26,36,27,37,28,38,29,48,19,100,...
    81,82,71,72,61,62,51,52,41,42,31,32,21,22,1,12,2,23,91,13,4,14,15,5,...
    16,6,17,7,8,18,10,9]';


% from the spreadsheet data for Dmitri's 3/24/08 PMd array (1024-0217)
elec2array = @(elecVec)(flipud(permute(reshape(elecVec,[10,10,3]),[2,1,3]))); 
% 3=RGB


% color the pins
coloration = zeros(max(pin2elec),3,max(unitIDs(:,2))); % 3=RGB
for iUnit = 1:size(unitIDs,1)
    thisElectrode = pin2elec(unitIDs(iUnit,1));
    coloration(thisElectrode,:,unitIDs(iUnit,2)) = C(iUnit,:);
end
coloration = sum(coloration,3)./sum(coloration~=0,3);
coloration(isnan(coloration)) = 0;

% paint onto the array
coloration = elec2array(coloration);
figure(fignum); 
subplot(2,3,3); cla;
imagesc(coloration);

end
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function colorElectrodeArrayByGroup(idx,unitIDs,fignum)
% color the electrode array according to the cluster ids


% from the spreadsheet data for Dmitri's 3/24/08 PMd array (1024-0217)
pin2elec = [93,92,94,95,75,96,85,97,86,98,87,88,77,99,66,89,76,90,67,79,...
    58,80,78,70,68,60,69,50,59,40,49,30,83,84,73,74,63,64,53,54,43,55,...
    44,45,33,46,34,65,24,56,35,47,25,57,26,36,27,37,28,38,29,48,19,100,...
    81,82,71,72,61,62,51,52,41,42,31,32,21,22,1,12,2,23,91,13,4,14,15,5,...
    16,6,17,7,8,18,10,9]';


% from the spreadsheet data for Dmitri's 3/24/08 PMd array (1024-0217)
elec2array = @(elecVec)(flipud(reshape(elecVec,[10,10])'));


% color the pins
coloration = zeros(max(pin2elec),max(unitIDs(:,2)));
for iUnit = 1:size(unitIDs,1)
    thisElectrode = pin2elec(unitIDs(iUnit,1));
    coloration(thisElectrode,unitIDs(iUnit,2)) = idx(iUnit);
end
coloration = sum(coloration,2)./sum(coloration~=0,2);
coloration(isnan(coloration)) = 0;

% paint onto the array
coloration = elec2array(coloration);
figure(fignum); 
subplot(2,3,2); cla;
imagesc(coloration)
% colorbar;

end
%-------------------------------------------------------------------------%


























