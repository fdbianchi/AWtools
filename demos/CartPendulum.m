
% AWtools: Cart-spring-pendulum example (plant LTI 2x1)
%
% From: G. Grimm, J. Hatfield, I. Postlethwaite, A. R.
%   Teel, M. C. Turner, and L. Zaccarian, “Antiwindup for stable linear
%   systems with input saturation: An LMI-based synthesis,” IEEE Trans.
%   Automat. Contr., vol. 48, no. 9, pp. 1509–1525, Sep. 2003.)
% 

% fbianchi - 2015-05-26
% fbianchi - 2021-07-22 - rev

% cleaning
clearvars
clc
close all

% =========================================================================
% Plant

A = [  0,      1,      0,    0; 
    -330.46, -12.15   -2.44  0;
       0,      0,      0,    1; 
    -812.61, -29.87, -30.10, 0];
Bu = [0; 2.71762; 0; 6.68268];
Bw = [0; 0; 0; 15.61];
C  = [1 0 0 0; 
      0 0 1 0];
Du = [0; 0];
Dw = [0; 0];
G = ss(A,[Bw Bu],C,[Dw, Du]);
G.y = 'y';
G.u = {'w','u'};

% LQG controller
K = [64.81, 213.12, 1242.27, 85.82];
L = [64, 2054,  -8, -1432;
     -8, -280, 142, 10169]';
Klqg = ss(A-Bu*K-L*C,L,K,0); 
Klqg.y = 'u';
Klqg.u = 'y';

% response test
% Gcl = connect(G,-Klqg,'w','y');

% =========================================================================
% AW design with the four options + simulations

sinFile = 'CartPendulum_sim.slx';

% Full order based on sectors 
opts.MinDamping = 0.1;
opts.MaxFreq = 100;
[Kaw1,Taw1,g1] = awsyn(G(:,2),Klqg,'fsec',[],opts);
Kaw = Kaw1; Taw = Taw1;
sim(sinFile)
t1  = yout.time;
y11 = yout.signals(1).values(:,3);
y12 = yout.signals(2).values(:,3);
u1  = yout.signals(3).values(:,3);

% Full order based on small gain theorm 
[Kaw2,Taw2,g2] = awsyn(G(:,2),Klqg,'fsgt',[],opts);
Kaw = Kaw2; Taw = Taw2;
sim(sinFile)
t2  = yout.time;
y21 = yout.signals(1).values(:,3);
y22 = yout.signals(2).values(:,3);
u2  = yout.signals(3).values(:,3);

% Static based on Tunner's paper
opts.MinBndU = 0.04;
opts.weigthY = 1e-6;
opts.weigthU = 1e-4;
[Kaw3,Taw3,g3] = awsyn(G(:,2),Klqg,'ssec',[],opts);
Kaw = Kaw3; Taw = Taw3;
sim(sinFile)
t3  = yout.time;
y31 = yout.signals(1).values(:,3);
y32 = yout.signals(2).values(:,3);
u3  = yout.signals(3).values(:,3);

% Low-order based on Tunner's paper
opts.MinBndU = 0.01;
opts.weigthY = 1e-6;
opts.weigthU = 1e-4;
F1 = tf(1,[1/200 1]);
F2 = ss(eye(2));
[Kaw4,Taw4,g4] = awsyn(G(:,2),Klqg,'psec',[],opts,F1,F2);
Kaw = Kaw4; Taw = Taw4;
sim(sinFile)
t4  = yout.time;
y41 = yout.signals(1).values(:,3);
y42 = yout.signals(2).values(:,3);
u4  = yout.signals(3).values(:,3);

% without saturation
y01 = yout.signals(1).values(:,1);
y02 = yout.signals(2).values(:,1);
u0  = yout.signals(3).values(:,1);

% with saturation
y01s = yout.signals(1).values(:,2);
y02s = yout.signals(2).values(:,2);
u0s  = yout.signals(3).values(:,2);

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

subplot(311); hold on
plot(t,y01, 'Color',0.7*[1 1 1],'linewidth',2);
plot(t,y01s,'Color',clines(1,:));
plot(t1,y11,'Color',clines(2,:));
plot(t2,y21,'Color',clines(3,:));
plot(t3,y31,'Color',clines(4,:));
plot(t4,y41,'Color',clines(5,:));
plot([t(1) t(end)],[0 0],'k--')
title('Comparison of different AW design methods in example of Cart-pendulum')
ylabel('p')
legend('without saturation','with saturation',...
    'AW fsec','AW fsgt','AW ssec','AW psec')

subplot(312); hold on
plot(t,y02, 'Color',0.7*[1 1 1],'linewidth',2);
plot(t,y02s,'Color',clines(1,:));
plot(t1,y12,'Color',clines(2,:));
plot(t2,y22,'Color',clines(3,:));
plot(t3,y32,'Color',clines(4,:));
plot(t4,y42,'Color',clines(5,:));
plot([t(1) t(end)],[0 0],'k--')
ylabel('\theta')

subplot(313); hold on
plot(t,u0, 'Color',0.7*[1 1 1],'linewidth',2);
plot(t,u0s,'Color',clines(1,:));
plot(t1,u1,'Color',clines(2,:));
plot(t2,u2,'Color',clines(3,:));
plot(t3,u3,'Color',clines(4,:));
plot(t4,u4,'Color',clines(5,:));
plot([t(1) t(end)],[0 0],'k--')
ylabel('u')
xlabel('time (s)')
