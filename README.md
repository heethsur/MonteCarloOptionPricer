# MonteCarloOptionPricer
A MonteCarlo Options Pricing Implementation in C#.

Takes option inputs and returns options price and associated greeks. 

Current feature set:
- Pricing European options; puts and calls
- Variance reduction technqiues for price paths including antithetic sampling and delta control variate
- Quasi-random using van der Corput for lower discrepancy and variance in price paths


Following developments are currently work in progress:
- Multithread execution for improved speed in generating price paths
- Alternative option payoffs 
