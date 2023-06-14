
% AWtools: Roll-Yaw Missile autopilot example (plant LTI of 2x2)
%
% From: "A new perspective on static and low order anti-windup synthesis", 
%   Turner & Postlethwaite, Int. Journal of Control, vol. 77, no. 1, 
%   pp. 27-44, 2004

% fbianchi - 2015-03-19
% fbianchi - 2021-07-22 - rev

% cleaning
clearvars
clc
close all

% =========================================================================
% Plant
Ap = [-0.818, -0.99900,  0.349; 
      80.290, -0.57900,  0.009;
   -2734.000,  0.05621, -2.100];
Bp = [0.147,     0.012;
   -194.400,    37.610;
  -2716.000, -1093.000];
Cp = [1 0 0; 0 1 0];
Dp = zeros(2);
P = ss(Ap,Bp,Cp,Dp);

% Controller
Ac1 = [-0.29, -107.80,   6.67,   -2.58,  -0.40; 
      107.68,  -97.81,  63.95,   -4.52,  -5.35;
       -6.72,   64.82, -54.19,  -40.79,   5.11; 
        3.21,    2.10,  29.56, -631.15, 429.89;
        0.36,   -3.39,   3.09, -460.03,  -0.74];
Bc1 = [2.28,   0.48;
     -40.75,   2.13;
      18.47,  -0.22;
      -2.07, -44.68;
      -0.98,  -1.18];
Cc1 = [0.86,  8.54,  -1.71, 43.91, 1.12;
       2.17, 39.91, -18.39, -8.51, 1.03];
Ac = [Ac1, Bc1; zeros(2,7)];
Bc = [zeros(5,2);eye(2)];
Cc = [Cc1, zeros(2)];
Dc = zeros(2);
C = ss(Ac,Bc,Cc,Dc);

% =========================================================================
% AW design with the four options + simulations

simfile = 'RollYawMissile_sim.slx';

% Full order based on sectors 
opts.MinDamping = 0.5;
opts.MaxFreq = 100;
[Kaw1,Taw1,g1] = awsyn(P,C,'fsec',[],opts);
Kaw = Kaw1; Taw = Taw1;
sim(simfile)
t1  = yout.time;
y11 = yout.signals(1).values(:,3);
y12 = yout.signals(2).values(:,3);
u11 = yout.signals(3).values(:,3);
u12 = yout.signals(4).values(:,3);

% Full order based on small gain theorm 
[Kaw2,Taw2,g2] = awsyn(P,C,'fsgt');%,[],opts);
Kaw = Kaw2; Taw = Taw2;
sim(simfile)
t2  = yout.time;
y21 = yout.signals(1).values(:,3);
y22 = yout.signals(2).values(:,3);
u21 = yout.signals(3).values(:,3);
u22 = yout.signals(4).values(:,3);

% Static based on Tunner's paper
opts.MinBndU = 1;
opts.weigthY = 1e-6;
[Kaw3,Taw3,g3] = awsyn(P,C,'ssec',[],opts);
Kaw = Kaw3; Taw = Taw3;
sim(simfile)
t3  = yout.time;
y31 = yout.signals(1).values(:,3);
y32 = yout.signals(2).values(:,3);
u31 = yout.signals(3).values(:,3);
u32 = yout.signals(4).values(:,3);

% Low-order based on Tunner's paper
opts.MinBndU = 0.03;
opts.weigthY = 1e-6;
opts.weigthU = 1e-2;
F1 = append(tf(2,[1 2]),1);
F2 = ss(eye(2));
[Kaw4,Taw4,g4] = awsyn(P,C,'psec',[],opts,F1,F2);
Kaw = Kaw4; Taw = Taw4;
sim(simfile)
t4  = yout.time;
y41 = yout.signals(1).values(:,3);
y42 = yout.signals(2).values(:,3);
u41 = yout.signals(3).values(:,3);
u42 = yout.signals(4).values(:,3);

% without sat
y01 = yout.signals(1).values(:,1);
y02 = yout.signals(2).values(:,1);
u01 = yout.signals(3).values(:,3);
u02 = yout.signals(4).values(:,3);
% with sat
y01s = yout.signals(1).values(:,2);
y02s = yout.signals(2).values(:,2);
u01s = yout.signals(3).values(:,3);
u02s = yout.signals(4).values(:,3);

t = yout.time;

fprintf('\n')
fprintf('--------------------------------------------------\n')
fprintf(' Performance:\n')
fprintf('\tFull order, sector boundeness: gamma = %7.3f\n',g1)
fprintf('\tFull order, small gain:        gamma = %7.3f\n',g2)
fprintf('\tStatic, sector boundeness:     gamma = %7.3f\n',g3)
fprintf('\tLow order, sector boundeness:  gamma = %7.3f\n',g4)
fprintf('--------------------------------------------------\n')
fprintf('\n')


%% ------------------------------------------------------------------------
% figures
figure('Position',[725    45   825   950]);
clines = lines(5);

subplot(411); hold on
plot(t,y01,'Color',0.7*[1 1 1],'linewidth',2); 
plot(t,y01s,'Color',clines(1,:));
plot(t1,y11,'Color',clines(2,:));
plot(t2,y21,'Color',clines(3,:));
plot(t3,y31,'Color',clines(4,:));
plot(t4,y41,'Color',clines(5,:));
plot([t(1) t(end)],[0 0],'k--')
title('Missile autopilot example: comparison of different AW methods')
ylabel('y_1')
legend('without saturation','with saturation',...
    'AW fsec','AW fsgt','AW ssec','AW psec')

subplot(412); hold on
plot(t,y02,'Color',0.7*[1 1 1],'linewidth',2); 
plot(t,y02s,'Color',clines(1,:));
plot(t1,y12,'Color',clines(2,:));
plot(t2,y22,'Color',clines(3,:));
plot(t3,y32,'Color',clines(4,:));
plot(t4,y42,'Color',clines(5,:));
plot([t(1) t(end)],[0 0],'k--')
ylabel('y_2')
axis([0 10 -15 10])

subplot(413); hold on
plot(t,u01,'Color',0.7*[1 1 1],'linewidth',2); 
plot(t,u01s,'Color',clines(1,:));
plot(t1,u11,'Color',clines(2,:));
plot(t2,u21,'Color',clines(3,:));
plot(t3,u31,'Color',clines(4,:));
plot(t4,u41,'Color',clines(5,:));
% plot([t(1) t(end)],-8*[1 1],'--','Color',0.7*[1 1 1])
% plot([t(1) t(end)], 8*[1 1],'--','Color',0.7*[1 1 1])
ylabel('u_1')

subplot(414); hold on
plot(t,u02,'Color',0.7*[1 1 1],'linewidth',2);
plot(t,u02s,'Color',clines(1,:));
plot(t1,u12,'Color',clines(2,:));
plot(t2,u22,'Color',clines(3,:));
plot(t3,u32,'Color',clines(4,:));
plot(t4,u42,'Color',clines(5,:));
plot([t(1) t(end)],-8*[1 1],'--','Color',0.7*[1 1 1])
plot([t(1) t(end)], 8*[1 1],'--','Color',0.7*[1 1 1])
ylabel('u_2')

xlabel('time (s)')
