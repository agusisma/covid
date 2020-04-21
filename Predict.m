%% Research code by Agus Hasan

clear;
clc;

%%
load JKT.txt; % load data: month | date | suspected | active cases | cummilative recovered | cummulative death
    
DATA = JKT;

%% Input for prediction
tp      = 30;  % prediction horizon
contact = 0.4; % contacts increase

%% Simulation
Ti = 9;
%%
tf  = length(DATA);
N   = sum(DATA(1,3:end));                    % number of population
CFR = DATA(end,end)/(sum(DATA(end,4:6)));    % case fatality rate
td  = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf);
tdp = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf+tp);
R0  = 3.2;
dt  = 0.01;
t   = dt:dt:tf;

%% Data matrix

C = [1 0 0 0 0;
     0 1 0 0 0; 
     0 0 1 0 0;
     0 0 0 1 0];

%% Initialization
xhat     = [N-1; 1; 0; 0; 0]; % initial condition
Pplus    = 0*eye(5);

%% Paramater
gamma  = (1-CFR)*(1/Ti);
kappa  = CFR*1/Ti;
alpha1 = 0.3;
sigma1 = 1.96; %95 CI
std_R  = 0.2;

%% Noise
QF = diag([10 10 10 10 std_R]);
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

xArray     = [];
xhatArray  = [];
    
for i=1:((tf-1)/dt)
     xhatArray = [xhatArray xhat]; 
     
     % prediction
     
     xhat(1) = xhat(1)-(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N;
     xhat(2) = xhat(2)+(gamma+kappa)*xhat(5)*xhat(1)*xhat(2)*dt/N-(gamma+kappa)*xhat(2)*dt;
     xhat(3) = xhat(3)+gamma*xhat(2)*dt;
     xhat(4) = xhat(4)+kappa*xhat(2)*dt;
     xhat(5) = xhat(5);

     xhat = xhat + sqrt(QF)*[randn randn randn randn randn]'*dt;
     
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

    y = y + sqrt(RF)*[randn randn randn randn]'*dt;
     
    Pmin  = FX*Pplus*FX'+QF;

    KF    = Pmin*C'*inv(C*Pmin*C'+RF);
%    QF    = alpha1*QF + (1-alpha1)*(KF*(y(:,i)-C*xhat)*(y(:,i)-C*xhat)'*KF');    

    % update 
    xhat  = xhat + KF*(y(:,i)-C*xhat);
    Pplus = (eye(5)-KF*C)*Pmin;
    
%    RF    = alpha1*RF + (1-alpha1)*((y(:,i)-C*xhat)*(y(:,i)-C*xhat)'+C*Pmin*C');
    xhat(5) = max(0,xhat(5)); % the reproduction number cannot be negative
end

%% Plotting

xhatArray(5,:) = filter(b,a,xhatArray(5,:));

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

%% Prediction
%% Initialization
xp   = [DATA(end,3); DATA(end,4); DATA(end,5); DATA(end,6)]; % initial condition

xpArray     = [];

mRt = mean(xhatRtArray(end-20:end));
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

xSpredic = [xhatSArray xpSArray];
xIpredic = [xhatIArray xpIArray];
xRpredic = [xhatRArray xpRArray];
xDpredic = [xhatDArray xpDArray];

%% Plotting

figure(1)
plot(tdp,xIpredic,'k','LineWidth',6)
%hold on
set(gca,'FontSize',24)
xline(datetime(2020,DATA(end,2),DATA(end,1)),'b','LineWidth',6)
text(tf-20,200,'\leftarrow Past','FontSize',24)
text(tf+10,200,'Future \rightarrow','FontSize',24)
grid on
grid minor

xhatRtArray(1,end)/R0
mean(xhatRtArray(end-7:end))/R0