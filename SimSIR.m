clear;
clc;

% by: Agus Hasan
% methods: https://www.pnas.org/content/pnas/115/50/12680.full.pdf
% https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf

% parameters
R0    = 3.0;       % basic reproduction number (https://www.ncbi.nlm.nih.gov/pubmed/32048815)

% Epidemic/pandemic is when R0>1. Current real-time estimate for R0 is between 2.8 and
% 3.3. Lockdown and social distancing could reduce R0 , hence flatten the
% infection curve

N     = 264000000; % total population (Indonesia)
Ni0   = 514;       % initial number of cases 22.03 (https://www.worldometers.info/coronavirus/country/indonesia/)
ti    = 10;        % mean of infection time (https://annals.org/aim/fullarticle/2762808/incubation-period-coronavirus-disease-2019-covid-19-from-publicly-reported)
B     = 0.06;      % critical beds requirement (https://www.ncbi.nlm.nih.gov/pubmed/32009128)
Ba    = 1.21;      % number of beds per 1000 of population (https://investor.id/national/rasio-bed-dibanding-populasi-di-indonesia-masih-rendah)
M     = 0.01;      % mortality rate (https://www.theguardian.com/world/2020/mar/22/what-is-coronavirus-and-what-is-the-mortality-rate)

% variance used in Monte Carlo simulation
varT  = 2;         % time incubation variance
varN  = 100;       % initial number of cases variance
varR  = 0.3;       % R0 variance

% for plotting
n     = 10000;     % number of simulation
nbins = 100;       % number plot for histogram
tf    = 300;       % simulation time (day)
t1    = 1:1:tf;
td    = datetime(2020,3,22) + caldays(1:tf);

y1A  = zeros(n,tf);
y2A  = zeros(n,tf);
y3A  = zeros(n,tf);
yPDF = zeros(n,1);

for i= 1:n
    
gamma = 1/(ti+varT*randn);
beta  = (R0+varR*randn)*gamma;
Ni    = Ni0+varN*randn;

tspan = [0 tf];
y0    = [1; Ni/N; 0];
[t,y] = ode45(@(t,y) odefcn(t,y,beta,gamma), tspan, y0);

yy1 = interp1(t,y(:,1),t1);
yy2 = interp1(t,y(:,2),t1);
yy3 = interp1(t,y(:,3),t1);

y1A(i,:) = yy1;
y2A(i,:) = yy2;
y3A(i,:) = yy3;

yPDF(i) = yy3(1,tf);

end

y1M = mean(y1A);
y2M = mean(y2A);
y3M = mean(y3A);

figure(1)
subplot(2,1,1)
plot(td,y1M,'b','LineWidth',3);
hold on
plot(td,y2M,'r','LineWidth',3);
hold on
plot(td,y3M,'g','LineWidth',3);
title('SIR Model');
xlabel('Date');
ylabel('Population fraction');
set(gca,'FontSize',24)
legend('Suspectible','Infectious','Recovered')
grid on;
grid minor;
subplot(2,1,2)
histogram(100*yPDF,nbins)
title('pdf of number of infected population');
xlabel('Percentage of infected population');
ylabel('Frequancy');
set(gca,'FontSize',24)
grid on
grid minor

figure(2)
subplot(2,1,1)
plot(td,B*N*y2M(1,1:tf)*(100000/N),'r','LineWidth',3);
hold on
plot(td,100*Ba*ones(tf,1),'g','LineWidth',3);
title('Number of available beds / 10^5 of population');
xlabel('Date');
set(gca,'FontSize',24)
legend('Needed','Available')
grid on;
grid minor;
subplot(2,1,2)
plot(td,[zeros(1,10) M*N*y2M(1,1:tf-10)*(100000/N)/ti],'b','LineWidth',3);
title('Mortality rate / 10^5 of population');
xlabel('Date');
set(gca,'FontSize',24)
grid on;
grid minor;