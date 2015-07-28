function [time,state,torque,evec] = simurobot(tfinal,ic,gains,constep,...
    contype,conparam,trajflag)
%SIMUROBOT: Simulates a controlled two link robot
%
% [t,y,tor,error] = simurobot(tf,ic,gains,constep,contype,conparam,trajflag)
%
% t - time vector
% y - state vector
% tor - torque vector
% error - controller errors at each update interval
% tf - final time for simulation in seconds
% ic - initial conditions = [theta 1;theta2;theta 1 dot;theta 2 dot]
% gains = [Kp,Kv,Ki], a 2x6 matrix where Kp,Kv,Ki are 2x2 gain matrices
% constep - controller update step size(seconds)
% contype - controller type: 1-Computed Torque, 2-PID, 3-Augmented PD
% conparam = [l1;l2;m1;m2;grav] - Controller Parameters
% trajflag = Trajectory Flag: 1-step, 2-sinusoid
%                                            Jeff Wendlandt, 10/12/94
t1 = clock;

statei = [ic(1);ic(2);ic(3);ic(4)];
state = [];
time = [];
torque = [];
evec = [];
dtime = [];
torint = [0;0];  % integrator storage
thetades = [];

param = conparam(1:5);
kp = gains(1:2,1:2);
kv = gains(1:2,3:4);
ki = gains(1:2,5:6);

% set up sinusoid trajectory
omega1 = 2*pi*0.5;
omega2 = 2*pi*1.0;
amp1 = 0.5;
amp2 = 0.4;

iold = 0;
i = 0
for i = 0:constep:(tfinal-constep),

    if ((i - iold) >= 0.2)
        iold = i;
        i
    end

    if trajflag == 1
        if i <= 0.2
            th1d = 0.0;
            th2d = 0.0;
            dth1d = 0.0;
            dth2d = 0.0;
            ddth1d = 0.0;
            ddth2d = 0.0;
        else
            th1d = 0.5;
            th2d = 0.9;
            dth1d = 0.0;
            dth2d = 0.0;
            ddth1d = 0.0;
            ddth2d = 0.0;
        end

    else
        if i <= 0.2
            th1d = 0.0;
            th2d = 0.0;
            dth1d = 0.0;
            dth2d = 0.0;
            ddth1d = 0.0;
            ddth2d = 0.0;
        else
            th1d = amp1*sin(omega1*(i-0.2));
            th2d = amp2*sin(omega2*(i-0.2));
            dth1d = amp1*omega1*cos(omega1*(i-0.2));
            dth2d = amp2*omega2*cos(omega2*(i-0.2));
            ddth1d = -amp1*omega1*omega1*sin(omega1*(i-0.2));
            ddth2d = -amp2*omega2*omega2*sin(omega2*(i-0.2));
        end

    end
    thetades = [thetades;th1d,th2d];

    th = statei(1:2);
    dth = statei(3:4);
    thd = [th1d;th2d];
    dthd = [dth1d;dth2d];
    ddthd = [ddth1d;ddth2d];

    e = th - thd;
    edot = dth - dthd;
    evec = [evec;e',edot'];
    dtime = [dtime;i];

    if contype == 1
        % Computed Torque
        massc = mass(th,param);
        corvec = coriotd(th,dth,param);
        potentfor = potfor(th,param);
        torint = -ki*e + torint;
        torquei = massc*(ddthd - kv*edot- kp*e + torint) + ...
            corvec + potentfor;

    elseif contype == 2
        % PID
        torint = -ki*e + torint;
        torquei = -kv*edot - kp*e + torint;

    else
        % Augmented PD
        massc = mass(th,param);
        potentfor = potfor(th,param);
        torquei = massc*ddthd + (corio(th,thd,param))*dthd + ...
            potentfor - kv*edot - kp*e;

    end

    % ode23 decreases cpu time in half
    %[temptime,tempstate] = ode45('robotsim',i,i+constep,[statei;torquei]);
    %[temptime,tempstate] = ode23('robotsim',i,i+constep,[statei;torquei]);
    [temptime,tempstate] = ode23('robotsim',[i i+constep],[statei;torquei]);

    len = size(temptime);
    len = len(1);

    time = [time;temptime(1:(len-1))];
    state = [state;tempstate(1:(len-1),1:4)];
    torque = [torque;tempstate(1:(len-1),5:6)];

    statei = tempstate(len,1:4)';


end
t2 = clock;
disp('Simulation Time is...');
secs = etime(t2,t1)

figure;
subplot(3,2,1),plot(time,state(:,1:2));
title('Theta 1 and 2 vs Time');
xlabel('Time(s)');
ylabel('Angle (radians)');

subplot(3,2,2),plot(time,torque);
title('Torque 1 and 2 vs Time');
xlabel('Time(s)');
ylabel('Torque (N-M)');

subplot(3,2,3),plot(dtime,evec(:,1:2));
title('Error 1 and 2 vs Time');
xlabel('Time(s)');
ylabel('Error (radians)');

subplot(3,2,4),plot(state(:,1),state(:,2));
title('Theta 2 vs Theta 1');
xlabel('Theta 1 (radians)');
ylabel('Theta 2 (radians)');
hold on;
plot(thetades(:,1),thetades(:,2),'b--');
hold off;

l1 = 4.0;
l2 = 3.0;
x = l1*cos(state(:,1)) + l2*cos(state(:,1) + state(:,2));
y = l1*sin(state(:,1)) + l2*sin(state(:,1) + state(:,2));
subplot(3,2,5), plot(x,y);
title('Trajectory of End Effector');
xlabel('X (m)');
ylabel('Y (m)');

subplot(3,2,6);
if contype == 1
    textt = 'Computed Torque';
elseif contype == 2
    textt = 'PID';
else
    textt = 'Augmented PD';
end
title(textt);
text(0.1,0.95,'gains are');
[g1,tmp] = sprintf('%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f \n',(gains(1,:)'));
[g2,tmp] = sprintf('%2.2f %2.2f %2.2f %2.2f %2.2f %2.2f \n',(gains(2,:)'));
text(0.1,0.8,g1);
text(0.1,0.7,g2);
text(0.1,0.5,'controller parameters are');
[g1,tmp] = sprintf('%2.2f %2.2f %2.2f %2.2f %2.2f',conparam);
text(0.1,0.4,g1);
text(0.1,0.2,'step size is');
text(0.1,0.1,num2str(constep));


