# covid

# SimSIR.m
This is an implemetation of SIR model with Monte Carlo simulation to predict the number of infected population due to covid-19.

R0 is called the basic reproduction number. Epidemic/pandemic is when R0>1. Current real-time estimate for R0 is between 2.8 and 3.3. Lockdown and social distancing could reduce R0, hence flatten the infection curve.

We assume that R0 has variance of 0.3 and the mean infection time is 10 days with variance of 2.

Figure 1 top shows the (normalize) number of suspectible, infectious, and recovered people. Figure 1 bottom shows the probability density function of infected population. Figure 2 top shows the number of critical beds care beds per 100000 of population. This can be used as an indication for the health care system capacity. Figure 2 bottom shows the mortality rate per 100000 of population.

# EKFRt
This program implements Extended Kalman Filter (EKF) to estimate the daily reproduction number (Rt). DATA.txt containts should containt the following information
month | date | suspected | active cases | cummilative recovered | cummulative death
