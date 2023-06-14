
% AWtools: F8 aircraft example (plant LTI 2x2)
%
% From: F. Wu and M. Soto, ‘Extended anti-windup control schemes for LTI 
%   and LFT systems with actuator saturations’, Int. J. of Robust and 
%   Nonlinear Control, vol. 14, no. 15, pp. 1255–1281, 2004

% fbianchi - 2015-03-19
% fbianchi - 2021-07-22 - rev

% cleaning
clearvars
clc
close all

% =========================================================================
% Plant
Ap = [-0.8, -0.0006, -12.00,   0;
       0,   -0.0140, -16.64, -32.2;
       1,   -0.0001,  -1.5,    0;
       1,    0,        0,      0];
Bp = [-19.00, -3;
       -0.66, -0.5;
       -0.16, -0.5;
        0,     0];
Cp = [0 0 0 1; 0 0 -1 1];
Dp = zeros(2);
G = ss(Ap,Bp,Cp,Dp,'InputName','u','OutputName','y');
umax = 15;

% -------------------------------------------------------------------------
% controller design
%
% weights
We = eye(2)*tf([0.5 2.5],[1 0.05]);
Wu = eye(2)*tf([1 5],[1 200]);
Wr = eye(2)*tf(4,[1 4]);
Wr.u = 'r'; Wr.y = 'yr';
%
% augmented plant
sbe = sumblk('e = r - y',2);
sbr = sumblk('v = y - yr',2);
Gau = connect(G,sbe,sbr,Wr,{'r','u'},{'v','u','e'});

% augmented plant + weights
Gauw = append(We,Wu,eye(2))*Gau;

% controller design
[K,CL,GAM,INFO] = hinfsyn(Gauw,2,2,'METHOD','lmi','DISPLAY','on');

% =========================================================================
% AW design with the four options + simulations

sinFile = 'F8_Example_sim.slx';

% Full order based on sectors 
opts.MinDecay = 0.03;
opts.MinDamping = 0.6;
opts.MaxFreq = 100;
[Kaw1,Taw1,g1] = awsyn(G,K,'fsec',[],opts);
Kaw = Kaw1; Taw = Taw1;
sim(sinFile)
t1  = yout.time;
y11 = yout.signals(1).values(:,3);
y12 = yout.signals(2).values(:,3);
u11 = yout.signals(3).values(:,3);
u12 = yout.signals(4).values(:,3);

% Full order based on small gain theorm 
opts.MinDecay = 0.01;
opts.MinDamping = 0.6;
opts.MaxFreq = 100;
[Kaw2,Taw2,g2] = awsyn(G,K,'fsgt',[],opts);
Kaw = Kaw2; Taw = Taw2;
sim(sinFile)
t2  = yout.time;
y21 = yout.signals(1).values(:,3);
y22 = yout.signals(2).values(:,3);
u21 = yout.signals(3).values(:,3);
u22 = yout.signals(4).values(:,3);

% Static based on Tunner's paper
opts.MinBndU = 1;
opts.weigthY = 1e-5;
[Kaw3,Taw3,g3] = awsyn(G,K,'ssec',[],opts);
Kaw = Kaw3; Taw = Taw3;
sim(sinFile)
t3  = yout.time;
y31 = yout.signals(1).values(:,3);
y32 = yout.signals(2).values(:,3);
u31 = yout.signals(3).values(:,3);
u32 = yout.signals(4).values(:,3);

% Low-order based on Tunner's paper
opts.MinBndU = 1;
opts.weigthY = 1e-4;
opts.weigthU = 1e-2;
F1 = append(tf(1,[1/10 1]),1);
F2 = ss(eye(2));
[Kaw4,Taw4,g4] = awsyn(G,K,'psec',[],opts,F1,F2);
Kaw = Kaw4; Taw = Taw4;
sim(sinFile)
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
fprintf('\tLow-order, sector boundeness:  gamma = %7.3f\n',g4)
fprintf('--------------------------------------------------\n')
fprintf('\n')

%% ------------------------------------------------------------------------
% figures
figure('Position',[725    45   825   950]);
clines = lines(5);

subplot(411); hold on
plot(t,y01, 'Color',0.7*[1 1 1],'linewidth',2);
plot(t,y01s,'Color',clines(1,:));
plot(t1,y11,'Color',clines(2,:));
plot(t2,y21,'Color',clines(3,:));
plot(t3,y31,'Color',clines(4,:));
plot(t4,y41,'Color',clines(5,:));
plot([t(1) t(end)],[0 0],'k--')
title('F8 example: comparison of different AW methods')
ylabel('y_1')
legend('without saturation','with saturation',...
    'AW fsec','AW fsgt','AW ssec','AW psec')

subplot(412); hold on
plot(t,y02, 'Color',0.7*[1 1 1],'linewidth',2);
plot(t,y02s,'Color',clines(1,:));
plot(t1,y12,'Color',clines(2,:));
plot(t2,y22,'Color',clines(3,:));
plot(t3,y32,'Color',clines(4,:));
plot(t4,y42,'Color',clines(5,:));
plot([t(1) t(end)],[0 0],'k--')
ylabel('y_2')

subplot(413); hold on
plot(t,u01, 'Color',0.7*[1 1 1],'linewidth',2);
plot(t,u01s,'Color',clines(1,:));
plot(t1,u11,'Color',clines(2,:));
plot(t2,u21,'Color',clines(3,:));
plot(t3,u31,'Color',clines(4,:));
plot(t4,u41,'Color',clines(5,:));
plot([t(1) t(end)],[0 0],'k--')
ylabel('u_1')

subplot(414); hold on
plot(t,u02, 'Color',0.7*[1 1 1],'linewidth',2);
plot(t,u02s,'Color',clines(1,:));
plot(t1,u12,'Color',clines(2,:));
plot(t2,u22,'Color',clines(3,:));
plot(t3,u32,'Color',clines(4,:));
plot(t4,u42,'Color',clines(5,:));
plot([t(1) t(end)],[0 0],'k--')
ylabel('u_2')

xlabel('time (s)')


