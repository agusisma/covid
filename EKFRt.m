%% Research code by Agus Hasan
% This code is used to estimate the value of daily reproduction number Rt
% based on Extended Kalman Filter (EKF) and low pass filter

clear;
clc;

%%
load DATA.txt; % load data: month | date | suspected | active cases | cummilative recovered | cummulative death

%%
tf  = length(DATA);
N   = sum(DATA(1,3:end));                    % number of population
CFR = DATA(end,end)/(sum(DATA(end,4:6)));    % case fatality rate
td  = datetime(2020,DATA(1,1),DATA(1,2)-1) + caldays(1:tf);
Ti  = 10;                                    % infection time

dt  = 0.01;
t   = dt:dt:tf;

%% Data matrix

C = [1 0 0 0 0;
     0 1 0 0 0; 
     0 0 1 0 0;
     0 0 0 1 0];

%% Initialization
xhat     = [N-1; 1; 0; 0; 0]; % initial condition 14.02
Pplus    = 0*eye(5);

%% Paramater
gamma  = (1-CFR)*(1/Ti);
kappa  = CFR*1/Ti;
sigma1 = 1.96; %95 CI
std_R  = 0.2;

%% Noise
QF = 1*eye(5);
RF = [100 0 0 0;0 10 0 0;0 0 1 0;0 0 0 std_R];

%% For plotting

windowSize = 500; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;

xArray     = [];
xhatArray  = [];

%% Simulation
for i=1:((tf-1)/dt)
     xhatArray = [xhatArray xhat]; 
     
     % prediction using discrete-time stochastic augmented compartmental model
     
     xhat(1) = xhat(1)-(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N;
     xhat(2) = xhat(2)+(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt;
     xhat(3) = xhat(3)+gamma*xhat(2)*dt;
     xhat(4) = xhat(4)+kappa*xhat(2)*dt;
     xhat(5) = xhat(5);

    % Extended Kalman filter
    % Calculating the Jacobian matrix
    FX    = [1-(gamma+kappa)*xhat(5)*xhat(2)*dt/N -(gamma+kappa)*xhat(5)*xhat(1)*dt/N 0 0 -(gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             (gamma+kappa)*xhat(5)*xhat(2)*dt/N 1+(gamma+kappa)*xhat(5)*xhat(1)*dt/N-(gamma+kappa)*dt 0 0 (gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             0 gamma*dt 1 0 0;
             0 kappa*dt 0 1 0;
             0 0 0 0 1];
    y = [interp1(0:1:tf-1,DATA(:,3),t,'makima');
         interp1(0:1:tf-1,DATA(:,4),t,'makima');
         interp1(0:1:tf-1,DATA(:,5),t,'makima');
         interp1(0:1:tf-1,DATA(:,6),t,'makima')];
    
    Pmin  = FX*Pplus*FX'+QF;

    KF    = Pmin*C'*inv(C*Pmin*C'+RF);
    
    % update 
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
ylabel('Active Cases')
set(gca,'FontSize',24)
legend('Estimation','Reported Cases')
grid on
grid minor
subplot(3,1,2)
plot(td,[xhatHArray DATA(end,5)],'LineWidth',6)
hold on
plot(td,DATA(:,5),'*r','LineWidth',6)
ylabel('Recovered')
set(gca,'FontSize',24)
grid on
grid minor
subplot(3,1,3)
plot(td,[xhatDArray DATA(end,6)],'LineWidth',6)
hold on
plot(td,DATA(:,6)','*r','LineWidth',6)
ylabel('Death')
xlabel('Date');
set(gca,'FontSize',24)
grid on
grid minor

H       = [xhatRArray xhatR];
curve1  = H + sigma1*std_R;
curve2  = max(H - sigma1*std_R,0);
x2      = [td, fliplr(td)];

figure(2)
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k');
alpha(0.5)
hold on;
plot(td,H,'k','LineWidth',6)
hold on
plot(td,ones(1,tf),'r','LineWidth',6)
title('Daily Reproduction Number (Rt)')
xlabel('Date');
set(gca,'FontSize',24)
grid on
grid minor
