% ee125 Lab 4
%-------------------------------------------------------------------------%
% Revised: 06/09/10
%   -cleaned up
% Created: Fall 2003
%   by JGM
%-------------------------------------------------------------------------%


clear;clc;

% gains
Kp = 11*eye(2);
Kv = 6*eye(2);
Ki = 0*eye(2);
Gain = [Kp Kv Ki];

%Two-link robot w/link lengths 4m and 3m, and masses 2kg and 3.5kg, resp.
roboparameters = [4 3 2 3.5 9.81]';      
IC = [0 0 0 0]';                            % Joints start at rest in 
                                            %  zero position

%---------------------------CT, step response-----------------------------%
update = 0.05;                              % Update controller at 20Hz
t_final = 3.0;                              % Simulate for 3.0s
[time state torque error] =...
    simurobot(t_final,IC,Gain,update,1,roboparameters,1);

%-----------------------------CT, figure 8--------------------------------%
pause;
t_final = 5.0;
[time, state, torque, error] =...
    simurobot(t_final,IC,Gain,update,1,roboparameters,2);

%-------------------decrease controller update frequency------------------%
pause;
t_final = 3.0; update = 0.25;
[time, state, torque, error] =...
    simurobot(t_final,IC,Gain,update,1,roboparameters,1);

%----------------------imperfect modelling of masses----------------------%
pause;
update = 0.05; roboparameters = [4;3;2.4;3.9;9.81];
[time, state, torque, error] =...
    simurobot(t_final,IC,Gain,update,1,roboparameters,1);

%----------------imperfect modelling, plus integral control---------------%
pause;
Ki = 0.5*eye(2);
Gain = [Kp Kv Ki];
[time, state, torque, error] =...
    simurobot(t_final,IC,Gain,update,1,roboparameters,1);

%------------------------------PID, step----------------------------------%
pause;
Gain = [3200 0 400 60 10 0; 0 700 60 100 0 1]; 
update = 0.01;
[time, state, torque, error] =...
    simurobot(t_final,IC,Gain,update,2,roboparameters,1);

%------------------------------PID figure 8-------------------------------%
pause;
t_final = 5.0;
[time, state, torque, error] =...
    simurobot(t_final,IC,Gain,update,2,roboparameters,2);

%-------------------------PID, step, double Kp----------------------------%
pause;
t_final = 3.0; Gain = [6400 0 400 60 10 0; 0 1400 60 100 0 1];
[time, state, torque, error] =...
    simurobot(t_final,IC,Gain,update,2,roboparameters,1);

%-------------------------PID, step, double Kv----------------------------%
pause;
Gain = [3200 0 800 120 10 0; 0 700 120 200 0 1];
[time, state, torque, error] =...
    simurobot(t_final,IC,Gain,update,2,roboparameters,1);

%-----------------------PID, step, Ki = O matrix--------------------------%
pause;
Gain = [3200 0 400 60 0 0; 0 700 60 100 0 0];
[time, state, torque, error] =...
    simurobot(t_final,IC,Gain,update,2,roboparameters,1);

%-----------------PID, step, decrease controller update-------------------%
pause;
Gain = [3200 0 400 60 10 0; 0 700 60 100 0 1]; update = 0.02;
[time, state, torque, error] =...
    simurobot(t_final,IC,Gain,update,2,roboparameters,1);
