% init
clear; clc; close all;
RESTART = 0;
filepre = 'results/hierwts/wts2DBP'; % wts1DPB


for ITER = 1:1% 100
    
    close all; clc;
    
    % load
    if RESTART
        % do some stuff to train a binbin1D network
        DBN
        filename = [filepre,'L1v',num2str(ITER,'%03.f'),'.mat'];
        save(filename,'wts','params');
    else
        % load results/finalwts/Std050.mat
        load('results\finalwts\wtsStandard140613','wts','params');
    end
    wts0 = wts;
    params0 = params;
    params0.smpls = 15;
    clear wts params;
    
    % set the new params
    params.MODEL = 'HierL2';
    params.NS = 'Joint-Angle-Left';
    params.t = params0.numsUnits(end);
    params.smin = params0.thmin;
    params.smax = params0.thmax;
    params.Nmods = 1; 
    params.mods = {'Joint-Angle-Left'};
    params.typeUnits = {'BP','Bernoulli'};
    params.numsUnits = [params.t+params0.numsUnits(1)/2, params.t];
    
    % copy over some old params
    params.Ndims = params0.Ndims;
    params.g = params0.g;
    params.swing = params0.swing;
    params.N = params0.N;
    params.C = params0.C;
    params.thmin = params0.thmin;
    params.thmax = params0.thmax;
    params.L1 = params0.L1;
    params.L2 = params0.L2;
    params.gst0 = params0.gst0;
    params.w = params0.w;
    params.q = params0.q;
    params.Ncases = params0.Ncases;
    params.respLength = params0.respLength;
    params.margin = params0.margin;
    params.gridsize = params0.gridsize;
    params.granularity = params0.granularity;
    params.nettype = params0.nettype;
    params.DBNmaxepoch = params0.DBNmaxepoch;
    params.epsilonw =  params0.epsilonw;
    params.epsilonvb = params0.epsilonvb;
    params.epsilonhb = params0.epsilonhb;
    params.weightcost = params0.weightcost;
    params.initialmomentum = params0.initialmomentum;
    params.finalmomentum = params0.finalmomentum;
    params.counterMax = params0.counterMax;
    params.numtestbatches = params0.numtestbatches;
    params.numCDsteps = params0.numCDsteps;
    
    
    
    % prepare the rbm
    numsUnits = params.numsUnits;
    numRBMs = length(numsUnits)-1;
    datagenargs = {'prevweights',wts0,'prevparams',params0};
    Nbatches = 1000;
    Ncases = params.Ncases;
    Nvis = params.numsUnits(1);
    wts = cell(numRBMs*2,1);
    
    
    % loop
    for i_rbm = 1:numRBMs
        Nhid = numsUnits(i_rbm+1);
        fprintf(1,'Pretraining Layer %i w/RBM: %d-%d \n',i_rbm,Nvis,Nhid);
        restart = 1;
        
        % train
        tic; EFH; toc;
        
        % pack together weights for saving (hid => recog., vis => gener.)
        wts{i_rbm} = [vishid; hidbiases];
        wts{numRBMs*2-i_rbm+1} = [vishid'; visbiases'];
        % filename = ['EncoderWtsFile',num2str(numsUnits(end))];
        filename = 'EncoderWtsFile';
        save(filename,'i_rbm','numsUnits','wts','params','epoch');
        
        % for next time through
        Nvis = Nhid;
        % batchdata = batchposhidprobs;
        
    end
    
    % save
    filename = [filepre,'L2v',num2str(ITER,'%03.f'),'.mat'];
    save(filename,'wts','params');
    
end

