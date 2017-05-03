function plotArmTracking(j,S0,roboparams,dynamics,varargin)
% plot the arm moving, and being tracked

% params
[Ncases,Ndims,Nmods,Nbatches] = size(S0);
dt = dynamics.A(1,Ndims+1);
setColors;
q = length(varargin)

% arm parameters
L1 = roboparams.armlengths(1);
L2 = roboparams.armlengths(2);
gstSHOULDERELBOW = [eye(3) [0 L1 0]'; 0 0 0 1];
gstELBOWHAND = [eye(3) [0 L1+L2 0]'; 0 0 0 1];
w1 = roboparams.w(:,1); q1 = roboparams.q(:,1);
w2 = roboparams.w(:,2); q2 = roboparams.q(:,2);


% figure(55); clf; hold on;
figure(56); clf; hold on;
for dim = 1:Ndims
    subplotHandle(dim) = subplot(1+Ndims,1,1+dim);
    hold on;
    for jj=1:q
        scatterHandle(dim,jj) = scatter(NaN,NaN);
    end
end
for t = 1:Nbatches
    
    % compute elbow and hand positions
    screwSHOULDER = screw(w1,cross(q1,w1),S0(j,1,2,t));
    screwELBOW = screw(w2,cross(q2,w2),S0(j,2,2,t));
    gstELBOWBASE = screwSHOULDER*gstSHOULDERELBOW;
    gstHANDBASE = screwSHOULDER*screwELBOW*gstELBOWHAND;
    
    elbow(1,1) = gstELBOWBASE(1,end);
    elbow(2,1) = gstELBOWBASE(2,end);
    
    hand(1,1) = gstHANDBASE(1,end);
    hand(2,1) = gstHANDBASE(2,end);
    
    
    % plot
    subplot(Ndims+1,1,1);
    hold on;
    plot(linspace(0,elbow(1)),linspace(0,elbow(2)),'m')
    plot(linspace(elbow(1),hand(1)),linspace(elbow(2),hand(2)),'m')
    axis tight; axis equal; axis([-45 30 -15 45]);
    
    % get predicted hand positions
    for fltr = 1:q
        Xhat = varargin{fltr}.Xhat;
        clr = varargin{fltr}.color;
        poshat = FK2link(Xhat(j,1:2,t),roboparams,1)';
        scatter(poshat(1),poshat(2),[],clr);
    end
    hold off;
    

    % set(0,'CurrentFigure',figure(56))
    for dim = 1:Ndims
        
        for fltr = 1:q
            thisErrs = squeeze(varargin{fltr}.Xhat(j,dim,1:t) - S0(j,2,dim,1:t));
            %%%% 2 gives prop---change this hack
            set(scatterHandle(dim,fltr),'XData',1:t,'YData',thisErrs,...
                'Marker','.','CData',varargin{fltr}.color);
        end
        
        hold off;
    end
    set(subplotHandle,'XLim',[0 Nbatches]);
    
    pause(dt);
    subplot(Ndims+1,1,1);
    cla
    
    
end
hold off;



end