# MonteCarloOptionPricer
A MonteCarlo Options Pricing Implementation in C#.

Takes option inputs and returns options price and associated greeks. 

Current feature set:
- A web API with interactive UI to access the pricer.
- Prices European, Asian, Digital, Barrier, and Lookback options. (Put and Call)
- Offers choice of variance reduction techniques that include:
    - Antithetic sampling
    - Control Variates*
- Quasi-random low-discrepancy Van der Corput sequence for generating random variables with improved variance
- Offers choice of multithreading for faster execution

*Control variate for European options is a delta-based control variate. All other options use underlying asset price as control variate.

Following developments are currently work in progress:
- Improved control variate for option types other than European.
 
