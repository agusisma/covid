%% Research code by Agus Hasan
% This code is used to estimate the value of daily reproduction number Rt
% based on Extended Kalman Filter (EKF) and low pass filter

clear;
clc;

%% load data
load IT.txt; % load data: date | month | susceptible | active cases | cummilative recovered | cummulative death

DATA = IT;
%% Infectious time
Tinf = 9;

%%
tp  = 30;                                    % prediction time
tf  = length(DATA);                          % simulation time
N   = sum(DATA(1,3:end));                    % number of population
CFR = DATA(end,end)/(sum(DATA(end,4:6)));    % case fatality rate
td  = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf);
tdp = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf+tp);
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

%% Prediction
cont = 0.05;
mRt  = mean(xhatRtArray(end-20:end));
R0   = max(xhatRtArray);

for m = 1:5

contact = cont*m;

xp   = [DATA(end,3); DATA(end,4); DATA(end,5); DATA(end,6)]; % initial condition
xpArray     = [];

cRt = contact*R0;
RtBArray = [];
for j=1:tp/dt
    RtB = mRt-((mRt-cRt)/(tp/dt))*j;
    RtBArray = [RtBArray RtB];
end
Rt = RtBArray;
for i=1:tp/dt
     xpArray = [xpArray xp]; 
     
     xp(1) = xp(1)-(gamma+kappa)*Rt(i)*xp(1)*xp(2)*dt/N;
     xp(2) = xp(2)+(gamma+kappa)*Rt(i)*xp(1)*xp(2)*dt/N-(gamma+kappa)*xp(2)*dt;
     xp(3) = xp(3)+gamma*xp(2)*dt;
     xp(4) = xp(4)+kappa*xp(2)*dt;
end

xpSArray  = [];
xpS       = xpArray(1,tf);
xpIArray  = [];
xpI       = xpArray(2,tf);
xpRArray  = [];
xpR       = xpArray(3,tf);
xpDArray  = [];
xpD       = xpArray(4,tf);
for i=1:tp
    xpSArray  = [xpSArray xpS];
    xpS       = xpArray(1,100*i);
    xpIArray  = [xpIArray xpI];
    xpI       = xpArray(2,100*i);
    xpRArray  = [xpRArray xpR];
    xpR       = xpArray(3,100*i);
    xpDArray  = [xpDArray xpD];
    xpD       = xpArray(4,100*i);
end

xIpredic(m,:) = [xhatIArray xpIArray];

end

% Plotting
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
title('Rt')
grid on
grid minor

figure(3)
plot(tdp,xIpredic(5,:),':c','LineWidth',6)
hold on;
plot(tdp,xIpredic(4,:),':g','LineWidth',6)
hold on;
plot(tdp,xIpredic(3,:),':k','LineWidth',6)
hold on;
plot(tdp,xIpredic(2,:),':b','LineWidth',6)
hold on;
plot(tdp,xIpredic(1,:),':r','LineWidth',6)
hold on;
set(gca,'FontSize',48)
xline(datetime(2020,DATA(end,2),DATA(end,1)),'b','LineWidth',6)
text(tf-20,0.3*DATA(end,4),'\leftarrow Past','FontSize',48)
text(tf,0.3*DATA(end,4),'Future \rightarrow','FontSize',48)
ylim([0 2*DATA(end,4)])
title('Projection')
grid on
grid minor

figure(4)
XRT = [1-(xhatRt/R0) xhatRt(end)/R0];
explode=[1 0];
h = pie(XRT,explode);
set(findobj(h,'type','text'),'fontsize',48)
title('Current PD Index')
set(gca,'FontSize',48)

figure(5)
subplot(2,2,1)
bar(td,[DATA(:,4) DATA(:,5) DATA(:,6)],'stacked')
set(gca,'color','none','FontSize',24)
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
title('COVID-19 Cases')
legend('Active Cases','Recovered','Death','Location','northwest')
grid on
grid minor

subplot(2,2,2)
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
set(gca,'color','none','FontSize',24)
ylim([0 6])
xlim([datetime(2020,DATA(1,2),DATA(1,1)), datetime(2020,DATA(end,2),DATA(end,1))])
legend('Confidence Interval 95%')
title('Real-Time Reproduction Number')
grid on
grid minor

subplot(2,2,3)
XRT = [1-(xhatRt/R0) xhatRt(end)/R0];
explode=[1 0];
h = pie(XRT,explode);
set(findobj(h,'type','text'),'fontsize',24)
title('Current Physical Distancing Index (PDI)')
set(gca,'FontSize',24)

subplot(2,2,4)
plot(tdp,xIpredic(5,:),':c','LineWidth',6)
hold on;
plot(tdp,xIpredic(4,:),':g','LineWidth',6)
hold on;
plot(tdp,xIpredic(3,:),':k','LineWidth',6)
hold on;
plot(tdp,xIpredic(2,:),':b','LineWidth',6)
hold on;
plot(tdp,xIpredic(1,:),':r','LineWidth',6)
set(gca,'color','none','FontSize',24)
xline(datetime(2020,DATA(end,2),DATA(end,1)),'b','LineWidth',6)
text(tf-20,0.3*DATA(end,4),'\leftarrow Past','FontSize',24)
text(tf,0.3*DATA(end,4),'Future \rightarrow','FontSize',24)
legend_str = {'PDI = 5%','PDI = 10%','PDI = 15%','PDI = 20%','PDI = 25%','Present'};
legend(legend_str(1:5),'Location','northwest')

ylim([0 2*DATA(end,4)])
title('30-Day Forecast')
grid on
grid minor