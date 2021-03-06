%% Research code by Agus Hasan
% This code is used to estimate the value of daily reproduction number Rt
% based on Extended Kalman Filter (EKF) and low pass filter

clear;
clc;

%% load data
load DATA.txt; % load data: date | month | susceptible | active cases | cummilative recovered | cummulative death

%% Infectious time
Tinf = 9;

%%
tf  = length(DATA);
N   = sum(DATA(1,3:end));                    % number of population
CFR = DATA(end,end)/(sum(DATA(end,4:6)));    % case fatality rate
td  = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf);

dt  = 0.01;
t   = dt:dt:tf;

%% Data matrix
C = [1 0 0 0 0;
     0 1 0 0 0; 
     0 0 1 0 0;
     0 0 0 1 0];
%% Parameters
sigma  = 1.96; %95 CI

%% Noise
QF = diag([10 10 10 10 0.2]);
RF = diag([100 10 10 1]);

%% For plotting
% Low pass filter
windowSize = 500; 
b          = (1/windowSize)*ones(1,windowSize);
a          = 1;
% For figure 2
curve11 = 0*ones(1,tf);
curve22 = 1*ones(1,tf);
x2      = [td, fliplr(td)];

%% Simulation
for j = 1:3
% Infection time
Ti     = Tinf-sigma+(j-1)*sigma;                     % infection time with CI 95%
gamma  = (1-CFR)*(1/Ti);
kappa  = CFR*1/Ti;

%% Initialization
xhat     = [N-1; 1; 0; 0; 0]; % initial condition
Pplus    = 0*eye(5);
xhatEff  = 0;

xArray       = [];
xhatArray    = [];
xhatEffArray = [];
    
for i=1:((tf-1)/dt)
     xhatArray    = [xhatArray xhat]; 
     xhatEffArray = [xhatEffArray xhatEff];      
     % adding reported data
     y = [interp1(0:1:tf-1,DATA(:,3),t,'makima');
         interp1(0:1:tf-1,DATA(:,4),t,'makima');
         interp1(0:1:tf-1,DATA(:,5),t,'makima');
         interp1(0:1:tf-1,DATA(:,6),t,'makima')];
     % predict
     xhat(1) = xhat(1)-(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N;
     xhat(2) = xhat(2)+(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt;
     xhat(3) = xhat(3)+gamma*xhat(2)*dt;
     xhat(4) = xhat(4)+kappa*xhat(2)*dt;
     xhat(5) = xhat(5);
    % calculating the Jacobian matrix
    FX    = [1-(gamma+kappa)*xhat(5)*xhat(2)*dt/N -(gamma+kappa)*xhat(5)*xhat(1)*dt/N 0 0 -(gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             (gamma+kappa)*xhat(5)*xhat(2)*dt/N 1+(gamma+kappa)*xhat(5)*xhat(1)*dt/N-(gamma+kappa)*dt 0 0 (gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             0 gamma*dt 1 0 0;
             0 kappa*dt 0 1 0;
             0 0 0 0 1];     
    Pmin  = FX*Pplus*FX'+QF;
    % update 
    KF    = Pmin*C'*inv(C*Pmin*C'+RF);
    xhat  = xhat + KF*(y(:,i)-C*xhat);
    Pplus = (eye(5)-KF*C)*Pmin;
    xhat(5) = max(0,xhat(5)); % the reproduction number cannot be negative
    xhatEff = (xhat(1)/N)*xhat(5); % calculating the effective repsoduction number
end

%% Plotting

xhatArray(5,:) = filter(b,a,xhatEffArray);

xhatSArray  = [];
xhatS       = xhatArray(1,tf);
xhatIArray  = [];
xhatI       = xhatArray(2,tf);
xhatRArray  = [];
xhatR       = xhatArray(3,tf);
xhatDArray  = [];
xhatD       = xhatArray(4,tf);
xhatRtArray = [];
xhatRt      = xhatArray(5,tf);
for i=1:tf-1
    xhatSArray  = [xhatSArray xhatS];
    xhatS       = xhatArray(1,100*i);
    xhatIArray  = [xhatIArray xhatI];
    xhatI       = xhatArray(2,100*i);
    xhatRArray  = [xhatRArray xhatR];
    xhatR       = xhatArray(3,100*i);
    xhatDArray  = [xhatDArray xhatD];
    xhatD       = xhatArray(4,100*i);
    xhatRtArray = [xhatRtArray xhatRt];
    xhatRt      = xhatArray(5,100*i);
end

xhatSArray  = [xhatSArray xhatS];
xhatIArray  = [xhatIArray xhatI];
xhatRArray  = [xhatRArray xhatR];
xhatDArray  = [xhatDArray xhatD];
xhatRtArray = [xhatRtArray xhatRt];

M(j,:) = xhatRtArray;

end

curve1      = M(1,:);
xhatRtArray = M(2,:);
curve2      = M(3,:);

%% Plotting
figure(1)
subplot(3,1,1)
plot(td,xhatIArray,'LineWidth',6)
hold on
plot(td,DATA(:,4),'*r','LineWidth',6)
ylabel('Active Cases')
set(gca,'FontSize',24)
legend('Estimated Cases','Reported Cases','Location','northwest')
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
grid on
grid minor
subplot(3,1,2)
plot(td,xhatRArray,'LineWidth',6)
hold on
plot(td,DATA(:,5),'*r','LineWidth',6)
ylabel('Recovered')
set(gca,'FontSize',24)
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
grid on
grid minor
subplot(3,1,3)
plot(td,xhatDArray,'LineWidth',6)
hold on
plot(td,DATA(:,6)','*r','LineWidth',6)
ylabel('Death')
set(gca,'FontSize',24)
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
grid on
grid minor

figure(2)
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k');
alpha(0.5)
hold on;
plot(td,xhatRtArray,'k','LineWidth',6)
hold on
inBetween = [curve11, fliplr(curve22)];
fill(x2, inBetween, 'g');
alpha(0.1)
hold on;
plot(td,ones(1,tf),'r','LineWidth',6)
set(gca,'FontSize',48)
ylim([0 6])
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
legend('Confidence Interval 95%')
grid on
grid minor

% RMS

RMSS = 0;
RMSI = 0;
RMSH = 0;
RMSD = 0;

for j = 1:tf
    RMSS = RMSS + sqrt(((xhatSArray(j)-DATA(j,3))/max(1,DATA(j,3)))^2);
    RMSI = RMSI + sqrt(((xhatIArray(j)-DATA(j,4))/max(1,DATA(j,4)))^2);
    RMSH = RMSH + sqrt(((xhatRArray(j)-DATA(j,5))/max(1,DATA(j,5)))^2);
    RMSD = RMSD + sqrt(((xhatDArray(j)-DATA(j,6))/max(1,DATA(j,6)))^2);
end
RMS  = (RMSS+RMSI+RMSH+RMSD)/tf
