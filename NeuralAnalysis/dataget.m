% DATADIR = '/Volumes/sabes/Documents/TODO/PEOPLE/ShenHelen/TwoTuning/Data';
DATADIR = 'C:\#DATA\Dmitri_datafiles';


DATES_D = {'D080616' 'D080728' 'D080822' 'D080617' 'D080729' 'D080826'...
           'D080619' 'D080731' 'D080827' 'D080620' 'D080820' 'D081013'};

DATES_E = {'E080201' 'E080205' 'E080207' 'E080208' 'E080212' 'E080213'...
           'E080214' 'E080215' 'E080226' 'E080227' 'E080228' 'E080229'...
           'E080305' 'E080306' 'E080308'};


%DATES = DATES_D;
% DATES = [DATES_D DATES_E];
DATES = DATES_D;
%DATES = {DATES_D{:} DATES_E{:}};
 
% for i=[1 3 4], figure(i); clf; end

Align    = 'Variable';
Epochs = {'EyeOn','Go'};


%% load data
for t=1:length(DATES)
    DATE = DATES{t};
    disp(DATE)
    
    
    for ep=1:size(Epochs,1)
        Epoch = Epochs(ep,:);
        
        TD(t,ep) = dataExtractor(DATE, Epoch, Align, 1);
    end
     
end
