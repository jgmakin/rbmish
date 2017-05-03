

clear; clc; close all;
load samplethings


% 
plot([dets00622; dets0125; dets025; dets050; dets100; dets250]'); % ; dets150; dets200; dets250]')
legend('N/8','N/4','N/2','N','2N','5N'); % ,'3N','4N','5N')
hold on;

opterrdet = 4.281587e-009;
plot(opterrdet*ones(1,16),'k:')
hold off;


% prop in det: 1.362495e-008
% prop out det for 5N network: 4.960751e-009

