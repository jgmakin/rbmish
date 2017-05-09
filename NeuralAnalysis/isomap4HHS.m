%%
clear; clc;

sessionDates = ['0616';'0617';'0619';'0620';'0728';'0729';'0820';'0822';...
    '0826';'0827';'1013'];
datadir = 'C:\#code\HHS\extracteddata\';

for iSession = 1:size(sessionDates,1)
   
    thisSessionDate = sessionDates(iSession,:);
    tag = ['D08',thisSessionDate]; 
    load([datadir,'KFtuningdataHHS',tag]);
    
    isomap4PNSdata(St,UnitSpikesT,str2double(thisSessionDate))
   
end
