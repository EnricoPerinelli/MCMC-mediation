# Monte Carlo Method for Assessing Mediation (MCMAM) function

Revised Monte Carlo simulation function for the computation of mediation effects (Monte Carlo Method for Assessing Mediation; MCMAM).

This function is based on the Preacher and Selig (2008) [web tool](https://www.quantpsy.org/medmc/medmc.htm), with the following enhancements:

- (a) Accepts **standard errors directly** (no need to manually square them),
- (b) Allows **setting a random seed** for reproducibility,
- (c) Facilitates **input of the covariance** between path a and path b as reported in the TECH3 output of Mplus, including scientific notation (number before and exponent after "D"),
- (d) Enables **automation of multiple mediation analyses** (e.g., cross-cultural studies with many countries or multi-group comparisons) by running the procedure via R code instead of manually using the online web application.