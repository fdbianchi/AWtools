

% AWtools: basic LPV example (plant LPV 1x1)

% fbianchi - 2021-07-29

% cleaning
clearvars; 
clc; 
close all

% ========================================================================
% Modelling

% parameter set
vert = [0 1 0;
        0 0 1];
pv = pset.Gral(vert);

% system matrix          
A(:,:,1) = [-2.0 1.0; -1.0 -2.0];
A(:,:,2) = [-1.3 0.5; -0.5 -1.3];
A(:,:,3) = [-3.0 1.0;  1.0 -3.0];
% all these matrices are the same at eache vertices
B = [1.0; 1.0];
C = [10 10];
D = 0;

% LPV model
pdG = ppss(A,B,C,D,pv);
pdG.u = 'u';    pdG.y = 'y';

% ========================================================================
% Control design

% Augmented plant
sb    = sumblk('e = r - y');
pdGau = connect(pdG,sb,{'r','u'},{'e','u','e'});
%
% weigths
W1 = tf(10,[1 0.01]);
W2 = tf([0.04 0.1],[0.004 1]);
Wout = append(W1,W2,1);
% augmented plant + weigths
pdGaw = Wout*pdGau;

% constraints
const(1) = synConst.Gain(1,1:2);
const(2) = synConst.Poles('MaxFreq',1000);

% synthesis
[pdK,constOut] = lpvsyn(pdGaw,3,2,const);
glpv = constOut(1).bound;

% =========================================================================
% AW design with the four options + simulations

simfile = 'BasicPolLPVmodel_Sat_sim.slx';

% Full order based on sectors 
opts.MinDamping = 0.8;
opts.MaxFreq = 1e3;
opts.weigthU = 1e-2;
[pdKaw1,pdTaw1,g1] = awsyn(pdG,pdK,'fsec',[],opts);
pdKaw = pdKaw1;
sim(simfile)
t1 = yout.time;
y1 = yout.signals(2).values(:,4);
u1 = yout.signals(3).values(:,3);

% Full order based on small gain theorm 
[pdKaw2,pdTaw2,g2] = awsyn(pdG,pdK,'fsgt',[],opts);
pdKaw = pdKaw2;
sim(simfile)
t2 = yout.time;
y2 = yout.signals(2).values(:,4);
u2 = yout.signals(3).values(:,3);

% Static based on Tunner's paper
opts.MinBndU = 1;
opts.weigthY = 1e-6;
[pdKaw3,pdTaw3,g3] = awsyn(pdG,pdK,'ssec',[],opts);
pdKaw = pdKaw3;
sim(simfile)
t3 = yout.time;
y3 = yout.signals(2).values(:,4);
u3 = yout.signals(3).values(:,3);

% Low-order based on Tunner's paper
opts.MinBndU = 0.03;
opts.weigthY = 1e-6;
opts.weigthU = 1e-2;
F1 = tf(1,[1/100 1]);
F2 = 1;
[pdKaw4,pdTaw4,g4] = awsyn(pdG,pdK,'psec',[],opts,F1,F2);
pdKaw = pdKaw4;
sim(simfile)
t4 = yout.time;
y4 = yout.signals(2).values(:,4);
u4 = yout.signals(3).values(:,3);

% common
r  = yout.signals(2).values(:,1);
p1 = yout.signals(1).values(:,1);
p2 = yout.signals(1).values(:,2);

% without sat
y0 = yout.signals(2).values(:,2);
u0 = yout.signals(3).values(:,1);
% with sat
ys = yout.signals(2).values(:,3);
us = yout.signals(3).values(:,2);

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
clines = lines(5);

% parameter trajectory
plot(pv)
hold on
plot(p1,p2,'Color',clines(1,:),'LineWidth',1)
xlabel('p_1'); xlabel('p_2')

% reponses
figure('Position',[725    45   825   720]);

subplot(211); hold on
plot([t(1) t(end)],[0 0],'k--')
hl(1) = plot(t,r,'k');
hl(2) = plot(t,y0,'Color',0.7*[1 1 1],'linewidth',2);
hl(3) = plot(t,ys,'Color',clines(1,:));
hl(4) = plot(t1,y1,'Color',clines(2,:));
hl(5) = plot(t2,y2,'Color',clines(3,:));
hl(6) = plot(t3,y3,'Color',clines(4,:));
hl(7) = plot(t4,y4,'Color',clines(5,:));
title('Basic Polytopic LPV model example: comparison of different AW methods')
ylabel('y')
legend(hl,'Reference','without saturation','with saturation',...
    'AW fsec','AW fsgt','AW ssec','AW psec')

subplot(212); hold on
plot(t,u0,'Color',0.7*[1 1 1],'linewidth',2);
plot(t,us,'Color',clines(1,:));
plot(t1,u1,'Color',clines(2,:));
plot(t2,u2,'Color',clines(3,:));
plot(t3,u3,'Color',clines(4,:));
plot(t4,u4,'Color',clines(5,:));
plot([t(1) t(end)],-0.2*[1 1],'--','Color',0.7*[1 1 1])
plot([t(1) t(end)], 0.2*[1 1],'--','Color',0.7*[1 1 1])
ylabel('u')

xlabel('time (s)')












