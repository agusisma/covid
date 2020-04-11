%% Research code by Agus Hasan

% This code is used to estimate the value of R0 based on Extended Kalman Filter

clear;
clc;

%%
load TUNR0.txt;

%%
tf  = length(TUNR0);
N   = sum(TUNR0(1,:));                    % number of population
CFR = TUNR0(end,4)/(sum(TUNR0(end,2:4))); % case fatality rate
td  = datetime(2020,1,30) + caldays(1:tf);
Ti  = 14;                                 % infection time

dt  = 0.01;
t   = dt:dt:tf;

%% Measurement matrix

C = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0];

%% Noise
QF = 0.01*eye(5);
RF = 1*eye(4);

%% Initialization
xhat     = [N-1; 1; 0; 0; 1]; % initial condition 14.02
Pplus    = 1*eye(5);

%% Paramater
gamma = (1-CFR)*(1/Ti);
kappa = CFR*1/Ti;

%% For plotting

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
    % calculating the Jacobian matrix
    FX    = [1-(gamma+kappa)*xhat(5)*xhat(2)*dt/N -(gamma+kappa)*xhat(5)*xhat(1)*dt/N 0 0 -(gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             (gamma+kappa)*xhat(5)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt 1+(gamma+kappa)*xhat(5)*xhat(1)*dt/N-(gamma+kappa)*dt 0 0 (gamma+kappa)*xhat(1)*xhat(2)*dt/N;
             0 gamma*dt 1 0 0;
             0 kappa*dt 0 1 0;
             0 0 0 0 1];

    Pmin  = FX*Pplus*FX'+QF;
    
    % update
    
    KF    = Pmin*C'*inv(C*Pmin*C'+RF);
    
    y = [interp1(0:1:tf-1,TUNR0(:,1),t);
         interp1(0:1:tf-1,TUNR0(:,2),t);
         interp1(0:1:tf-1,TUNR0(:,3),t);
         interp1(0:1:tf-1,TUNR0(:,4),t)];
     
    xhat  = xhat + KF*(y(:,i)-C*xhat);
    Pplus = (eye(5)-KF*C)*Pmin;
end

%% Plotting

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
subplot(2,1,1)
plot(td,[xhatSArray TUNR0(end,1)],'LineWidth',6)
hold on
plot(td,TUNR0(:,1),'*r','LineWidth',6)
ylabel('S')
set(gca,'FontSize',24)
legend('Estimasi','Data Lapangan')
grid on
grid minor
subplot(2,1,2)
plot(td,[xhatIArray TUNR0(end,2)],'LineWidth',6)
hold on
plot(td,TUNR0(:,2),'*r','LineWidth',6)
ylabel('I')
set(gca,'FontSize',24)
grid on
grid minor

figure(2)
subplot(2,1,1)
plot(td,[xhatHArray TUNR0(end,3)],'LineWidth',6)
hold on
plot(td,TUNR0(:,3),'*r','LineWidth',6)
ylabel('R')
set(gca,'FontSize',24)
legend('Estimasi','Data Lapangan')
grid on
grid minor
subplot(2,1,2)
plot(td,[xhatDArray TUNR0(end,4)],'LineWidth',6)
hold on
plot(td,TUNR0(:,4)','*r','LineWidth',6)
ylabel('D')
xlabel('Tanggal');
set(gca,'FontSize',24)
grid on
grid minor

figure(3)
semilogy(td,[xhatRArray xhatRArray(end)],'LineWidth',6)
hold on
semilogy(td,ones(1,tf),'LineWidth',6)
title('Estimasi Nilai R0')
ylabel('R0');
xlabel('Tanggal');
set(gca,'FontSize',24)
grid on
grid minor
