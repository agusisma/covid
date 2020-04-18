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
Ti  = 9;                                    % infection time

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
alpha1 = 1;
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
     
     % prediction
     
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
%    QF    = alpha1*QF + (1-alpha1)*(KF*(y(:,i)-C*xhat)*(y(:,i)-C*xhat)'*KF');    

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

xhatSArray = [xhatSArray xhatS];
xhatIArray = [xhatIArray xhatI];
xhatHArray = [xhatHArray xhatH];
xhatDArray = [xhatDArray xhatD];
xhatRArray = [xhatRArray xhatR];

figure(1)
subplot(3,1,1)
plot(td,xhatIArray,'LineWidth',6)
hold on
plot(td,DATA(:,4),'*r','LineWidth',6)
ylabel('Active Cases')
set(gca,'FontSize',24)
legend('Estimation','Reported Cases')
grid on
grid minor
subplot(3,1,2)
plot(td,xhatHArray,'LineWidth',6)
hold on
plot(td,DATA(:,5),'*r','LineWidth',6)
ylabel('Recovered')
set(gca,'FontSize',24)
grid on
grid minor
subplot(3,1,3)
plot(td,xhatDArray,'LineWidth',6)
hold on
plot(td,DATA(:,6)','*r','LineWidth',6)
ylabel('Death')
xlabel('Date');
set(gca,'FontSize',24)
grid on
grid minor

curve1  = xhatRArray + sigma1*std_R;
curve2  = max(xhatRArray - sigma1*std_R,0);
x2      = [td, fliplr(td)];

figure(2)
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k');
alpha(0.5)
hold on;
plot(td,xhatRArray,'k','LineWidth',6)
hold on
plot(td,ones(1,tf),'r','LineWidth',6)
title('Daily Reproduction Number (Rt)')
xlabel('Date');
set(gca,'FontSize',24)
grid on
grid minor

% RMS

% RMSS = sqrt(((xhatSArray)'-DATA(1:end-1,3))'*((xhatSArray)'-DATA(1:end-1,3)))/length(xhatSArray);
% RMSI = sqrt(((xhatIArray)'-DATA(1:end-1,4))'*((xhatIArray)'-DATA(1:end-1,4)))/length(xhatIArray);
% RMSR = sqrt(((xhatHArray)'-DATA(1:end-1,5))'*((xhatHArray)'-DATA(1:end-1,5)))/length(xhatHArray);
% RMSD = sqrt(((xhatDArray)'-DATA(1:end-1,6))'*((xhatDArray)'-DATA(1:end-1,6)))/length(xhatDArray);
% RMS  = RMSS+RMSI+RMSR+RMSD

for j = 1:tf
    RMSS = sqrt(((xhatSArray(j)-DATA(j,3))/DATA(j,3))^2);
    RMSI = sqrt(((xhatIArray(j)-DATA(j,4))/DATA(j,4))^2);
    RMSH = sqrt(((xhatHArray(j)-DATA(j,5))/DATA(j,5))^2);
    RMSD = sqrt(((xhatDArray(j)-DATA(j,6))/DATA(j,6))^2);
    RMS  = (RMSS+RMSI+RMSH+RMSD)*100                        % Percentage of RMS
end
