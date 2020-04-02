%% Research code by Agus Hasan

% This code is to estimate the value of R0 based on Extended Kalman Filter
% The data were taken from 14.02 until 30.03

clear;
clc;

%%
load TUNR01.txt; % Assuming CFR 1%
load TUNR02.txt; % Assuming CFR 4%

%%
tf  = 45;
dt  = 0.01;
t   = dt:dt:tf;
td  = datetime(2020,2,14) + caldays(1:tf);
N   = 264000000;                % number of population
CFR = 0.04;                     % death rate 4%
Ti  = 10;                       % infection time

%% Measurement matrix

C = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0];

%% Noise
QF = 0.001*eye(5);
RF = 0.001*eye(4);

%% Initialization
xhat     = [N-1; 1; 0; 0; 1]; % initial condition 14.02
Pplus    = 0*eye(5);

%% Paramater
gamma = (1-CFR)*(1/Ti);
kappa = CFR*1/Ti;

%% For plotting

xArray     = [];
xhatArray  = [];

%%
% Simulation
for i=1:(tf/dt)
     xhatArray = [xhatArray xhat]; 
     
     % prediction
     
     xhat(1) = xhat(1)-(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N;
     xhat(2) = xhat(2)+(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt;
     xhat(3) = xhat(3)+gamma*xhat(2)*dt;
     xhat(4) = xhat(4)+kappa*xhat(2)*dt;
     xhat(5) = xhat(5);

    % Extended Kalman filter
    % calculating the Jacobian matrix
    FX    = [1-(gamma+kappa)*xhat(5)*xhat(2)*dt/N -(gamma+kappa)*xhat(5)*xhat(1)*dt/N 0 0 -(gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             (gamma+kappa)*xhat(5)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt 1+(gamma+kappa)*xhat(5)*xhat(1)*dt/N-(gamma+kappa)*dt 0 0 (gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             0 gamma*dt 1 0 0;
             0 kappa*dt 0 1 0;
             0 0 0 0 1];

    Pmin  = FX*Pplus*FX'+QF;
    
    % update
    
    KF    = Pmin*C'*inv(C*Pmin*C'+RF);
    
    % insert data/measurement
%     y = [interp1(0:1:tf,TUNR01(:,1),t);
%          interp1(0:1:tf,TUNR01(:,2),t);
%          interp1(0:1:tf,TUNR01(:,3),t);
%          interp1(0:1:tf,TUNR01(:,4),t)];
    y = [interp1(0:1:tf,TUNR02(:,1),t);
         interp1(0:1:tf,TUNR02(:,2),t);
         interp1(0:1:tf,TUNR02(:,3),t);
         interp1(0:1:tf,TUNR02(:,4),t)];
     
    xhat      = xhat + KF*(y(:,i)-C*xhat);
    Pplus = (eye(5)-KF*C)*Pmin;
end

figure(1)
subplot(4,1,1)
plot(t,xhatArray(1,:),'LineWidth',3)
hold on
%plot(0:1:tf,TUNR01(:,1),'*r','LineWidth',3)
plot(0:1:tf,TUNR02(:,1),'*r','LineWidth',3)
ylabel('S')
title('Model SIRD')
set(gca,'FontSize',24)
grid on
grid minor
subplot(4,1,2)
plot(t,xhatArray(2,:),'LineWidth',3)
hold on
%plot(0:1:tf,TUNR01(:,2),'*r','LineWidth',3)
plot(0:1:tf,TUNR02(:,2),'*r','LineWidth',3)
ylabel('I')
set(gca,'FontSize',24)
grid on
grid minor
subplot(4,1,3)
plot(t,xhatArray(3,:),'LineWidth',3)
hold on
%plot(0:1:tf,TUNR01(:,3),'*r','LineWidth',3)
plot(0:1:tf,TUNR02(:,3),'*r','LineWidth',3)
ylabel('R')
set(gca,'FontSize',24)
grid on
grid minor
subplot(4,1,4)
plot(t,xhatArray(4,:),'LineWidth',3)
hold on
%plot(0:1:tf,TUNR01(:,4),'*r','LineWidth',3)
plot(0:1:tf,TUNR02(:,4),'*r','LineWidth',3)
ylabel('D')
xlabel('Date');
set(gca,'FontSize',24)
grid on
grid minor

figure(2)
plot(t,xhatArray(5,:),'LineWidth',3)
hold on
plot(t,mean(xhatArray(5,2000:end))*ones(1,4500),'r','LineWidth',3)
title('Estimasi Nilai R0')
ylabel('R0');
xlabel('Date');
legend('R0 dari EKF','Rata2 setelah steady state')
set(gca,'FontSize',24)
grid on
grid minor
