function RFshifts(wts,params)
% This function plots the tuning curves (for position) of the coordinate-
% transformation model under two different gaze angles, as well as
% hypothetical retinotopic tuning curves at the second gaze angle, given
% the tuning curves at the first gaze angle.
%
% This is supposed to correspond to the figure in your paper, the code for
% which somehow got lost--hence the need to recreate it at this late date.
%
% Presently, it just prints tuning curves from random hidden units, so some
% assemby is still required.  Some day you may change it to select
% partially shifting, fully shifting, non-shifting, and anti-shifting
% tuning curves.
%
% NB: hard-coded for length(params.mods) = 3, params.Ndims = 1;

%-------------------------------------------------------------------------%
% Created: 06/30/14
%   by JGM
%-------------------------------------------------------------------------%


% params
if checkGPUavailability, dataclass = 'gpuArray'; else dataclass = 'double'; end
xmin = params.roboparams.posmin;
xmax = params.roboparams.posmax;
emin = params.roboparams.eyemin;
emax = params.roboparams.eyemax;
Ndims = params.Ndims;
Nmods = length(params.mods);
hidDstrbs = params.typeUnits{2};
hidNums = params.numsUnits{2};

Nexamples = 1000;
fractionsOfGazeSpace = [1/4 3/4];
gazeAngles = scalefxn(fractionsOfGazeSpace,zeros(2,1),ones(2,1),...
    emin*ones(2,1),emax*ones(2,1));
clrs = ['r','g','b'];
%%%% thresh = 0.15;


% prepare figure
Nrows = 4; Ncols = 4;
figure(13); p = panel(); p.pack(Nrows,Ncols);
% inds = (1:Nrows*Ncols);
inds = [1,85,63,87,99,16,28,90,73,74,140,79,27,113,20,92];
    

% [~,inds] = sort(rand(Nhids,1));
% inds = inds(1:Nrows*Ncols);

% stimuli (mostly placeholders)
S = NaN(Nexamples,Ndims,Nmods);
S(1:Nexamples,1:Ndims,strcmp(params.mods,'Hand-Position')) =...
    linspace(xmin,xmax,Nexamples)';

% for each eye position
for iGaze = 1:length(fractionsOfGazeSpace)
    
    % fixed gaze angle, all the different x's, and the corresponding joints
    e = gazeAngles(iGaze)*ones(Nexamples,1);
    x = S(:,:,strcmp(params.mods,'Hand-Position'));
    th = IK2link(x - e,params.roboparams,1);
    
    % store
    S(:,:,strcmp(params.mods,'Joint-Angle')) = th;
    S(:,:,strcmp(params.mods,'Gaze-Angle')) = e;
    
    % compute the multisensory tuning curves
    params.typeUnits{1} = 'Dirac';
    gains = mean([params.gmin; params.gmax]);
    params.gmin = gains;
    params.gmax = gains;
    R = generateData(Nexamples,params,'stimuli',S);
    V = invParamMap(R,wts{1}(1:end-1,:),wts{1}(end,:),hidDstrbs,hidNums,params);
    
    
    
    % if  max(v) > thresh
    % end
    for j = 1:Ncols
        for i = 1:Nrows
            p(i,j).select();
            hold on;
            theseX = S(:,:,strcmp(params.mods,'Hand-Position')); 
            plot(theseX,V(:,inds(i+(j-1)*Nrows)),'color',clrs(iGaze),...
                'Linewidth',2.0);
            
            if iGaze==1
                % V1 is the response to (X-E) at E=e1; call it phi(x).  For
                % retinotopic neurons, the response at E=e2 is the response
                % at X-e1+e1-e2 = X-e1 - (e2-e1), i.e. phi(X-(e2-e1)) where
                % the parenthetical term is the gaze shift. So the response
                % shifts *rightward* for a positive gaze shift (your case).
                % (If the old peak was at X=10 and the shift is 5, the new
                % peak will be at X=15.)  So we *add* the gaze shift to the
                % x axis.
                plot(theseX+diff(gazeAngles),V(:,inds(i+(j-1)*Nrows)),...
                    'color',clrs(iGaze),'Linestyle',':','Linewidth',2.0);
                axis off
                
            end

            axis([emin+xmin emax+xmax 0 1]);
            if i~=Nrows, set(gca, 'xtick', []); end
            if j~=1, set(gca, 'ytick', []); end
            hold off;
        end
    end
    


end
p.de.margin = 6;



end