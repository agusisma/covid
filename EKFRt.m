%% Research code by Agus Hasan
% This code is used to estimate the value of daily reproduction number Rt
% based on Extended Kalman Filter (EKF)

clear;
clc;

%%
load DATA.txt; % load data: month | date | suspected | active cases | cummilative recovered | cummulative death

%%
tf  = length(DATA);
N   = sum(DATA(1,3:end));                    % number of population
CFR = DATA(end,end)/(sum(DATA(end,4:6)));    % case fatality rate
td  = datetime(2020,DATA(1,1),DATA(1,2)) + caldays(1:tf);
Ti  = 14;                                    % infection time

dt  = 0.01;
t   = dt:dt:tf;

%% Data matrix

C = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0];

%% Noise
QF = 1*eye(5);
RF = [10 0 0 0;0 10 0 0;0 0 5 0; 0 0 0 2];

%% Initialization
xhat     = [N-1; 1; 0; 0; 1]; % initial condition 14.02
Pplus    = 1*eye(5);

%% Paramater
gamma = (1-CFR)*(1/Ti);
kappa = CFR*1/Ti;

%% For plotting

windowSize = 500; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

xArray     = [];
xhatArray  = [];

%% Simulation
for i=1:((tf-1)/dt)
     xhatArray = [xhatArray xhat]; 
     
     % prediction
     
     xhat(1) = xhat(1)-(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N;
     xhat(2) = xhat(2)+(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt;
     xhat(3) = xhat(3)+gamma*xhat(2)*dt;
     xhat(4) = xhat(4)+kappa*xhat(2)*dt;
     xhat(5) = xhat(5);

    % Extended Kalman filter
    % Calculating the Jacobian matrix
    FX    = [1-(gamma+kappa)*xhat(5)*xhat(2)*dt/N -(gamma+kappa)*xhat(5)*xhat(1)*dt/N 0 0 -(gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             (gamma+kappa)*xhat(5)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt 1+(gamma+kappa)*xhat(5)*xhat(1)*dt/N-(gamma+kappa)*dt 0 0 (gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             0 gamma*dt 1 0 0;
             0 kappa*dt 0 1 0;
             0 0 0 0 1];

    Pmin  = FX*Pplus*FX'+QF;
    
    % update
    
    KF    = Pmin*C'*inv(C*Pmin*C'+RF);
    
    y = [interp1(0:1:tf-1,DATA(:,3),t);
         interp1(0:1:tf-1,DATA(:,4),t);
         interp1(0:1:tf-1,DATA(:,5),t);
         interp1(0:1:tf-1,DATA(:,6),t)];
     
    xhat  = xhat + KF*(y(:,i)-C*xhat);
    Pplus = (eye(5)-KF*C)*Pmin;
end

%% Plotting

xhatArray(5,:) = filter(b,a,xhatArray(5,:));

xhatSArray = [];
xhatS      = xhatArray(1,tf);
xhatIArray = [];
xhatI      = xhatArray(2,tf);
xhatHArray = [];
xhatH      = xhatArray(3,tf);
xhatDArray = [];
xhatD      = xhatArray(4,tf);
xhatRArray = [];
xhatR      = xhatArray(5,tf);
for i=1:tf-1
    xhatSArray = [xhatSArray xhatS];
    xhatS      = xhatArray(1,100*i);
    xhatIArray = [xhatIArray xhatI];
    xhatI      = xhatArray(2,100*i);
    xhatHArray = [xhatHArray xhatH];
    xhatH      = xhatArray(3,100*i);
    xhatDArray = [xhatDArray xhatD];
    xhatD      = xhatArray(4,100*i);
    xhatRArray = [xhatRArray xhatR];
    xhatR      = xhatArray(5,100*i);
end

figure(1)
subplot(3,1,1)
plot(td,[xhatIArray DATA(end,4)],'LineWidth',6)
hold on
plot(td,DATA(:,4),'*r','LineWidth',6)
ylabel('Kasus Aktif')
set(gca,'FontSize',24)
grid on
grid minor
subplot(3,1,2)
plot(td,[xhatHArray DATA(end,5)],'LineWidth',6)
hold on
plot(td,DATA(:,5),'*r','LineWidth',6)
ylabel('Sembuh')
set(gca,'FontSize',24)
legend('Estimasi','Data Lapangan')
grid on
grid minor
subplot(3,1,3)
plot(td,[xhatDArray DATA(end,6)],'LineWidth',6)
hold on
plot(td,DATA(:,6)','*r','LineWidth',6)
ylabel('Meninggal')
xlabel('Tanggal');
set(gca,'FontSize',24)
grid on
grid minor

H= [xhatRArray xhatRArray(end)];
x = 1:numel(H);
sigma = 1.96; %95 CI

figure(2)
std_dev = 0.25;
curve1 = H + sigma*std_dev;
curve2 = H - sigma*std_dev;
x2 = [td, fliplr(td)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k');
alpha(0.5)
hold on;
plot(td,H,'k','LineWidth',6)
hold on
plot(td,ones(1,tf),'r','LineWidth',6)
title('Estimasi Angka Reproduksi Harian (Rt)')
xlabel('Tanggal');
set(gca,'FontSize',24)
grid on
grid minor
