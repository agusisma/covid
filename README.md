# covid
This is an implemetation of SIR model with Monte Carlo simulation to predict the number of infected population.

An interesting parameter is R0, which is the basic reproduction number. Epidemic/pandemic is when R0>1. Current real-time estimate for R0 is between 2.8 and 3.3. Lockdown and social distancing could reduce R0, hence flatten the infection curve.

We assume the R0 has variance of 0.3 and the mean infection time is 10 days with variance of 2.

Figure 1 top shows the (normalize) number of suspectible, infectious, and recovered people. Figure 1 bottom shows the probability density function of infected population. Figure 2 top shows the number of critical beds care beds per 100000 of population. This can be used asn an indication for the capability of health care system. Figure 2 bottom shows the mortality rate per 100000 of population.
