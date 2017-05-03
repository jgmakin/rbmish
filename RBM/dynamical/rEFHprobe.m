function [MSE,MSE0,MSE1] = rEFHprobe(wts,params,varargin)
% rEFHprobe  Probe how well a recurrent EFH has learned
%
% USAGE:
%   load([getdir('data'),'RBMish\EFHs\wts_2DrEFHbb_160304.mat'])
%   [MSE,MSE0,MSE1] = rEFHprobe(wts,params)
%
%   load([getdir('data'),'RBMish/EFHs/wts_MoCap_23-Nov-2016.mat']);
%   params.testortrain = 'test';
%   [MSE,MSE0,MSE1] = rEFHprobe(wts,params)
%
%   load('C:\#DATA\RBMish\EFHs\wts_polyphonicmusic_120_20-Dec-2016.mat')
%   params.testortrain = 'valid';
%   [MSE,MSE0,MSE1] = rEFHprobe(wts,params);
%
% rEFHprobe computes next-frame MSE, baseline MSEs, and confabulates
% trajectories from the recurrent EFH defined by wts.  It works for
% (R)TRBMs as well as rEFHs.


%-------------------------------------------------------------------------%
% Revised: 12/05/16
%   -added functions for MoCap data
% Revised: 01/20/16
%   -added function for writing animated gifs
% Created: 08/10/15
%   by JGM
%-------------------------------------------------------------------------%

if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end


% how should the confabulation be initialized?
GIBBSINITIALIZE = 0;

% test next frame errors, also checking 0th- and 1st-order models
if ~isempty(varargin)
    R = varargin{1};    X = varargin{2};    Q = varargin{3};
else
    [R,X,Q] = params.getTestData(dataclass);
end
[MSE0, MSE1] = getBaselineMSEs(shortdata(floor(size(R,1)/Q.T),3,R),params);
[MSE, Zbar] = testEFHNextFrameError(R,X,Q,wts,params);

% how should you initialize?
if GIBBSINITIALIZE, Zf = GibbsInit(wts,params); else Zf = Zbar(1,:,end); end

keyboard

% generate backwards
switch params.datatype
    case 'bouncingballs' % ?,'ECoG',
        T = 400;
        Y = backwardGenerate(T,Zf,wts{2}(1:end-1,:),wts{2}(end,:),0.1,0,params);
        dispBouncingBalls(Y,1,params);
        writeGifTraj(Y,1,['confab_',params.datatype,'_1'],params);
    case 'MoCap'
        T = 100;
        Y = backwardGenerate(T,Zf,wts{2}(1:end-1,:),wts{2}(end,:),0.1,0,params);
        [skel, Motion] = dispMoCap(Y,1);
        writeMoCapMovie(skel, Motion{1}, 1/30);
    case 'polyphonicmusic'
        T = 50;
        Y = backwardGenerate(T,Zf,wts{2}(1:end-1,:),wts{2}(end,:),0,0.9,params);
        [y,Fs] = playPolyphonicMusic(Y,1);
        
        y = 0.95*y/max(abs(y));
        audiowrite('rEFHmusic.wav', y, Fs);
    case 'HDFRE'
        
        T = 1000;
        Y = backwardGenerate(T,Zf,wts{2}(1:end-1,:),wts{2}(end,:),0.5,0,params);
        charset = params.charset;        
        switch params.typeUnits{1}{1}
            case 'Bernoulli'
                %%% this has some problems....
                %%%inds = longdata(round(Y))*(2.^[0:6]') + 1;
                inds = round(longdata(Y)*(2.^[0:6]')) + 1;
                inds(inds > length(charset)) = length(charset);
            case 'Categorical'
                inds = round(longdata(Y))*(1:(length(charset)-1))' + 1;
                %[~,iii] = max(longdata(Y),[],2);
                %zeroinds = sum(round(longdata(Y)),2)==0;
                %inds(zeroinds) = 0;
                %inds = inds+1;
        end
        charset(inds)
        
    case 'spikecounts'
        T = 1000;
        Y = backwardGenerate(T,Zf,wts{2}(1:end-1,:),wts{2}(end,:),0.5,0,params);
        
        [Y,Z,XhatStatic,XhatDynamic] = fakah(T,Zf,wts{2}(1:end-1,:),...
            wts{2}(end,:),0.2,0,wts,params);
        
        Y = longdata(Y);
        Z = longdata(Z);
        X = XhatStatic;
        
        keyboard
        
end

return




% generate in the forward direction
Yhat2 = forwardGenerate(Zbar,43,0,'confabhidstates',wts,params);
dispBouncingBalls(Yhat2,3,params)



end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
function Zf = GibbsInit(wts,params)

Ncdsteps = 50;

% wts
vishid = wts{1}(1:end-1,:);
hidbiases = wts{1}(end,:);
visbiases = wts{2}(end,:)';

% params
hidDstrbs   = params.typeUnits{2};
hidNums     = params.numsUnits{2};
visDstrbs   = [hidDstrbs, params.typeUnits{1}];
visNums     = [hidNums,   params.numsUnits{1}];

% initial hidden state
visInit = zeros(1,sum(visNums));
hidmeans = invParamMap(visInit,zeros(size(vishid),'like',vishid),...
    zeros(size(hidbiases),'like',vishid),hidDstrbs,hidNums,params);
hidInit = sampleT(hidmeans,hidDstrbs,hidNums,params);

% now Gibbs sample
[~, hidmeans] = CDstepper(hidInit,vishid,visbiases,...
    hidbiases,hidDstrbs,visDstrbs,hidNums,visNums,Ncdsteps,params);

%%% Here you're assume that CDstepper's final output is means! (this is
%%% commented into place for training recurrent EFH's....)
Zf = sampleT(hidmeans,hidDstrbs,hidNums,params);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function Y = backwardGenerate(T,Zf,hidvis,visbiases,noiseZZ,noiseYZ,params)

%%% It seems *theoretically* quite interesting that the hidbiases are not
%%% required for backwards generation---except getting Zf.

%%%
% this only works for two-layer networks....
%%%


% visible units and distributions
hidDstrbs   = params.typeUnits{2};
hidNums     = params.numsUnits{2};
visDstrbs   = [hidDstrbs, params.typeUnits{1}];
visNums     = [hidNums,   params.numsUnits{1}];
inputInd    = sum(hidNums)+1;



% malloc
Y = zeros(size(Zf,1),params.numsUnits{1},T,'like',Zf);

% loop backwards through time
for t = T:-1:1
    Vbar = invParamMap(Zf,hidvis,visbiases,visDstrbs,visNums,params);
    sampleDstrbs = visDstrbs;
    
    % noiseless transition? noiseless emission?
    if rand > noiseZZ, sampleDstrbs{1} = 'Dirac'; end
    if rand > noiseYZ, [sampleDstrbs{2:end}] = deal('Dirac'); end
    
    % now sample (or possibly "sample")
    V = sampleT(Vbar,sampleDstrbs,visNums,params);
    Y(:,:,t) = V(:,inputInd:end);
    Zf = V(:,1:(inputInd-1));
    
end
fprintf('done generating\n');

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function dispBouncingBalls(Y,iTraj,params)

for t = 1:size(Y,3)
    topogrids = displayshape(Y(iTraj,:,t),params);
    imagesc(topogrids{1})
    title(num2str(t))
    pause(0.02);
end

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function writeGifTraj(Y,iTraj,filename,params)

filename = [filename,'.gif'];

imagestack = ones(params.N+2,params.N+2,1,size(Y,3),'like',Y);
for t = 1:size(Y,3)
    topogrids = displayshape(Y(iTraj,:,t),params);
    imagestack(2:end-1,2:end-1,:,t) = topogrids{1};
end
imagestack = uint8((2^8-1)*gather(imagestack));
newmap = gray(256);
imwrite(imagestack,newmap,filename,'gif','LoopCount',Inf,'DelayTime',0.02);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [skel,Motion] = dispMoCap(Y,iTraj)

% function dispMoCap(Y,params)
%%%% not a great way to do this
%%%
load([getdir('data'),'MoCap/data'],'skel','Motion');
for ii=1:length(Motion),
    Motion{ii} = cast(Motion{ii}(1:4:end,:),'like',Y);
end
%%%
% load([getdir('data'),'MoCap/Jog1_M'],'skel','Walking');
% skel.type = 'mit';
% for ii=1:length(Walking),
%     Motion{ii} = cast(Walking{ii}(1:4:end,:),'like',Y);
% end
%%%

data = mocapBase2EFHdata(skel,Motion);
data.seqindex = mat2cell(1:size(Y,1)*size(Y,3),1,size(Y,3)*ones(size(Y,1),1));
data.R = longdata(Y);
%%%

Motion = mocapEFHdata2Base(skel,data);
figure(201); expPlayData(skel, Motion{iTraj}, 1/30);

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [MSE0, MSE1] = getBaselineMSEs(R,params)

if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end

% hard-coded
t0 = 1;

% 0th-order
MSE0 = mean(mean(longdata(diff(R,[],3)).^2));

% 1st-order
%%%%%% FIX ME!!
params.testortrain = 'train';
[R,X,Q] = params.getTestData(dataclass);
% I doubt this will work!  getTestData's dependence on params.testortrain
% was set when it was created!!  Rather, you need to just use
% params.getLatents and params.getData.
%%%%%%
R = shortdata(floor(size(R,1)/Q.T),3,R);
Rp = longdata(R(:,:,1:(end-1)));
Rf = longdata(R(:,:,2:end));
B = linregress(Rp,Rf);
SE = (longdata(R(:,:,t0:(end-1)))*B - longdata(R(:,:,(t0+1):end))).^2;
MSE1 = mean(SE(:));

end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function writeMoCapMovie(skel, channels, frameLength, xlim, ylim, zlim)
% modified very mildly from expPlayData.m
%   by JGM
% 
% Version 1.000 
%
% Code provided by Graham Taylor, Geoff Hinton and Sam Roweis 
%
% For more information, see:
%     http://www.cs.toronto.edu/~gwtaylor/publications/nips2006mhmublv
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with
% a note saying that the original programs are available from our
% web page.
% The programs and documents are distributed without any warranty, express or
% implied.  As the programs were written for research purposes only, they have
% not been tested to the degree that would be advisable in any important
% application.  All use of these programs is entirely at the user's own risk.
%
% Based on skelPlayData.m version 1.1
% Copyright (c) 2006 Neil D. Lawrence
%
% We support two types of skeletons:
%  1) Those built from the CMU database (acclaim)
%     http://mocap.cs.cmu.edu/
%  2) Those built from data from Eugene Hsu (mit)
%     http://people.csail.mit.edu/ehsu/work/sig05stf/
% EXPPLAYDATA Play skel motion capture data.
% Data is in exponential map representation
%
% Usage: [xlim, ylim, zlim] = expPlayData(skel, channels, frameLength)

if nargin < 3
    frameLength = 1/120;
end
clf

handle = expVisualise(channels(1, :), skel);

if nargin < 6
    %We didn't specify the limits of the motion
    %So calculate the limits

        xlim = get(gca, 'xlim');
        minY1 = xlim(1);
        maxY1 = xlim(2);
        ylim = get(gca, 'ylim');
        minY3 = ylim(1);
        maxY3 = ylim(2);
        zlim = get(gca, 'zlim');
        minY2 = zlim(1);
        maxY2 = zlim(2);
        for ii = 1:size(channels, 1)
            Y = exp2xyz(skel, channels(ii, :));
            minY1 = min([Y(:, 1); minY1]);
            minY2 = min([Y(:, 2); minY2]);
            minY3 = min([Y(:, 3); minY3]);
            maxY1 = max([Y(:, 1); maxY1]);
            maxY2 = max([Y(:, 2); maxY2]);
            maxY3 = max([Y(:, 3); maxY3]);
        end
        xlim = [minY1 maxY1];
        ylim = [minY3 maxY3];
        zlim = [minY2 maxY2];
end

set(gca, 'xlim', xlim, ...
    'ylim', ylim, ...
    'zlim', zlim);

% Play the motion
yrvid = VideoWriter('mocapvideo.mp4', 'MPEG-4');
yrvid.FrameRate = 1/frameLength;
open(yrvid)
for jj = 1:size(channels, 1)
    expModify(handle, channels(jj, :), skel);
    F = getframe;
    writeVideo(yrvid,F);
end
close(yrvid);


end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
function [y, Fs] = playPolyphonicMusic(Y,iSeq)

% init
Fs = 44100;
notelength = 0.3;   % s (that's what B-L has in his code)
Note0 = 21;         % you just have to know this
path(path,genpath(('C:\#code\miditools')));
%%% don't worry, it won't add it "twice"

pianoroll = squeeze(Y(iSeq,:,:));
notes = pianoRoll2matrix(pianoroll, notelength, Note0 + (0:(size(Y,2)-1)));
mididream = matrix2midi(notes);
y = midi2audio(mididream, Fs, 'sine');
soundsc(y,Fs);

end
%-------------------------------------------------------------------------%