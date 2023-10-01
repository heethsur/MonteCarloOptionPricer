using System;

namespace MonteCarloOptionPricing // MonteCarlo Pricing with variance reduction techniques
{
    class Program
    {
        static void Main(string[] args)
        {
            

 
            EuropeanOption option = new();
               
            Console.WriteLine("Welcome to MonteCarlo Option Pricer!");

            Console.WriteLine("Enter option type (C for Call, P for Put):");
            option.OptionType = Console.ReadLine()?.ToUpper();

            Console.WriteLine("Enter underlying price:");
            option.UnderlyingPrice = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter strike price:");
            option.StrikePrice = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter time to expiration (in years):");
            option.TimeToExpiration = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter volatility:");
            option.Volatility = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter risk-free rate:");
            option.RiskFreeRate = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter whether to use antithetic sampling (true for Yes, false for No):");
            option.useAntitheticSampling = Convert.ToBoolean(Console.ReadLine());

            Console.WriteLine("Enter whether to use delta-based control variate (true for Yes, false for No):");
            option.useControlVariate = Convert.ToBoolean(Console.ReadLine());

                    
            // Output of results
            Console.WriteLine("------");
            Console.WriteLine("European " + option.OptionType  + " Option with underlying price " + option.UnderlyingPrice + " and strike " + option.StrikePrice + " expiring in " + option.TimeToExpiration + " years"  );
            Console.WriteLine("------");
            Console.WriteLine("Option price: " + option.OptionPrice()[0]);
            Console.WriteLine("------");
            option.FetchGreeks();
            Console.WriteLine("------");
            Console.WriteLine("Standard error: " + option.OptionPrice()[1]);
        }
    }

    class EuropeanOption
    // This class contains all properties for a European Option, and methods for pricing and associated greeks
    {
        public string? OptionType { get; set; } // Call or Put?
        public double UnderlyingPrice { get; set; } // Spot Price
        public double StrikePrice { get; set; } // Strike Price
        public double TimeToExpiration { get; set; } // t2x
        public double Volatility { get; set; } // sigma
        public double RiskFreeRate { get; set; } // mu
        public bool useAntitheticSampling {get; set; } // using Antithetic Sampling?
        public bool useControlVariate {get; set; } // using Control Variate?
        private int TimeSteps = 252; // timesteps are fixed at 252
        private int Simulations = 100000; // number of simulations are fixed at 100,000

        public double [] OptionPrice() // this method produces Option Price and its Standard Error using Monte Carlo
        {
            double [] result = new double[2]; // this result array will hold Option Price, and Standard Error
            double dt = TimeToExpiration/TimeSteps; 
            int beta = -1; // take opposite to delta position in underlying
            double sum_Payoffs = 0; // will hold total payoffs 
            double sum_Payoffs_sq = 0; // will hold total payoffs squared

            for(int i = 1; i < (Simulations + 1); i++) 
            {// for each simulation generate a 252 period path of underlying prices following a GBM
                double S_t1 = UnderlyingPrice;
                double S_t2 = UnderlyingPrice * (useAntitheticSampling? 1:0); // will hold Antithetic complements
                double cv1 = 0; // control variate
                double cv2 = 0; // will hold Antithetic control variate
                double [] Z = Norm.RandomValues1(TimeSteps); // generate 252 random variables ~ N(0,1) along each path

                for(int j = 1; j < (TimeSteps + 1); j++) // for each path
                {
                    double t = (j-1)*dt; // time elapsed along the path
                    double d_1 = delta(S_t1,t); // delta during jth day in the path
                    double d_2 = delta(S_t2,t); // '' for anithetic complement on jth day
                    double S_t1_n = S_t1*Math.Exp((RiskFreeRate - (Volatility*Volatility/2))*dt + (Volatility*Math.Sqrt(dt)*Z[j-1]));
                    double S_t2_n = S_t2*Math.Exp((RiskFreeRate - (Volatility*Volatility/2))*dt + (Volatility*Math.Sqrt(dt)*-Z[j-1]))*(useAntitheticSampling? 1:0); // Antithetic complement
                    cv1 += d_1*(S_t1_n - S_t1 * Math.Exp(RiskFreeRate*dt))*(useControlVariate? 1:0); //delta controle variate
                    cv2 += d_2*(S_t2_n - S_t2 * Math.Exp(RiskFreeRate*dt))*(useControlVariate? 1:0)*(useAntitheticSampling? 1:0); //''f for antithetic
                    S_t1 = S_t1_n; // update t+1 price as start for t+2 price
                    S_t2 = S_t2_n; // '' for antithetic complement
                }
                double Payoff = 0;
                if (OptionType == "C") // Terminal Option Payoff (+ Hedge Cashflows if control variate)
                {// Payoff dependent on whether used Antithetic or no
                    Payoff = 0.5*(Math.Max(0, S_t1 - StrikePrice) + beta*cv1 + (useAntitheticSampling? (Math.Max(0, S_t2 - StrikePrice) + beta*cv2):(Math.Max(0, S_t1 - StrikePrice) + beta*cv1)));
                }
                else
                {
                    Payoff = 0.5*(Math.Max(0, StrikePrice - S_t1) + beta*cv1 + (useAntitheticSampling? (Math.Max(0, StrikePrice - S_t2) + beta*cv2):(Math.Max(0, StrikePrice - S_t1) + beta*cv1)));
                }
                
                sum_Payoffs += Payoff; 
                sum_Payoffs_sq += Payoff*Payoff;
            }

            result[0] = Math.Round(((sum_Payoffs/Simulations) * Math.Exp(-RiskFreeRate*TimeToExpiration)),3) ; // discount terminal payoffs
            double SD = Math.Sqrt((sum_Payoffs_sq - sum_Payoffs*sum_Payoffs/Simulations)*Math.Exp(-2*RiskFreeRate*TimeToExpiration)/(Simulations-1)); // Standard Deviation
            result[1] = SD/Math.Sqrt(Simulations); // Standard error
            return result;
        }
        
        //Greeks calculated similarly to Project 1 from last semester
        private double d1(double S, double t2x) // takes Underlying and Time to Expiration as inputs since delta needs to be iteratively calculated for each time step
        // compute d1
        {
            double v = (Math.Log(S/StrikePrice) + (RiskFreeRate + (Volatility*Volatility)/2)*t2x)/(Volatility * Math.Sqrt(t2x));
            return v;
        }

        private double d2(double S, double t2x)
        // compute d2
        {
            double v = d1(S,t2x) - Volatility*Math.Sqrt(t2x);
            return v;
        }

        private double delta(double S,double t2x)
        // compute option delta
        {
            double v = 0;
            if (OptionType == "C")
            {
                v = Norm.Cdf(d1(S,t2x));
            }
            else
            {
                v = Norm.Cdf(d1(S,t2x)) - 1;
            }
            return Math.Round(v,3);
        }

        private double gamma()
        // compute option gamma
        {
            double v = Norm.Pdf(d1(UnderlyingPrice,TimeToExpiration))/(UnderlyingPrice*Volatility*Math.Sqrt(TimeToExpiration));
            return Math.Round(v,3);
        }

        private double vega()
        // compute option vega
        {
            double v = UnderlyingPrice*Norm.Pdf(d1(UnderlyingPrice,TimeToExpiration))*Math.Sqrt(TimeToExpiration);
            return Math.Round(v,3);
        }

        private double theta()
        // compute option theta
        {
            double v = 0;
            if (OptionType == "C")
            {
                v = (-UnderlyingPrice*Norm.Pdf(d1(UnderlyingPrice,TimeToExpiration))*Volatility)/(2*Math.Sqrt(TimeToExpiration)) - RiskFreeRate*StrikePrice*Math.Exp(-RiskFreeRate*TimeToExpiration)*Norm.Cdf(d2(UnderlyingPrice,TimeToExpiration));
            }
            else
            {
                v = (-UnderlyingPrice*Norm.Pdf(d1(UnderlyingPrice,TimeToExpiration))*Volatility)/(2*Math.Sqrt(TimeToExpiration)) + RiskFreeRate*StrikePrice*Math.Exp(-RiskFreeRate*TimeToExpiration)*Norm.Cdf(-d2(UnderlyingPrice,TimeToExpiration));
            }
            return Math.Round(v,3);
        }

        private double rho()
        // compute option rho
        {
            double v = 0;
            if (OptionType == "C")
            {
                v = StrikePrice*TimeToExpiration*Math.Exp(-RiskFreeRate*TimeToExpiration)*Norm.Cdf(d2(UnderlyingPrice,TimeToExpiration));
            }
            else
            {
                v = -StrikePrice*TimeToExpiration*Math.Exp(-RiskFreeRate*TimeToExpiration)*Norm.Cdf(-d2(UnderlyingPrice,TimeToExpiration));
            }
            return Math.Round(v,3);
        }

        public void FetchGreeks() //Function to get all Greeks at once
        {
            Console.WriteLine("Option Greeks:");
            Console.WriteLine("------");
            Console.WriteLine("Delta: " + delta(UnderlyingPrice,TimeToExpiration));
            Console.WriteLine("Gamma: " + gamma());
            Console.WriteLine("Vega: " + vega());
            Console.WriteLine("Theta: " + theta());
            Console.WriteLine("Rho: " + rho());
            
        }

    }

    
    class Norm
    {
        // This class contains everything to do with Normal Distribution: compute pdf, cdf, and generate standard normal random variables
        public static double Pdf(double x)
        // function to calculate normal pdf of x
        {
            return Math.Exp(-0.5 * x * x) / Math.Sqrt(2 * 3.141592653589793238462643383279502884197);
        }

        public static double Cdf(double x)
        // approximation of cdf for x since closed form solution is not tractable. Refer to Marsagila polar method
        {
            double t = 1.0 / (1.0 + 0.2316419 * Math.Abs(x));
            double b1 = 0.319381530;
            double b2 = -0.356563782;
            double b3 = 1.781477937;
            double b4 = -1.821255978;
            double b5 = 1.330274429;
            double t2 = t * t;
            double t3 = t2 * t;
            double t4 = t3 * t;
            double t5 = t4 * t;
            double phi = Pdf(x);
            double cdf = 1.0 - phi * (b1 * t + b2 * t2 + b3 * t3 + b4 * t4 + b5 * t5);
            if (x < 0.0)
            {
                cdf = 1.0 - cdf;
            }
            return cdf;
        }

        
        private static double[] GenerateVanDerCorputSequence(int base_value)
        {// Generate VanDerCprput Sequence given a base value
            double[] sequence = new double[1000]; 

            for (int i = 0; i < 1000; i++)
            {
                int denom = base_value;
                double value = 0;

                int current = i+1;  

                while (current > 0)
                {
                    value += (current % base_value) / (double)denom;
                    current /= base_value;
                    denom *= base_value;
                }

                sequence[i] = value;
            }

            return sequence;
        }

        public static double[] RandomValues(int count)
        {  // Borrowed from Project 4 last semester: generating a standard normal variable using polar rejection
            Random random = new();
            double [] sample = GenerateVanDerCorputSequence(2);
            double[] values = new double[count];

            for (int i = 0; i < count; i++)
            {
                int randomIndex1 = random.Next(0, sample.Length);
                int randomIndex2 = random.Next(0, sample.Length);
                double u1 = sample[randomIndex1];
                double u2 = sample[randomIndex2];
                // transform uniform distribution interval to [-1,1]
                double v1 = 2 * u1 - 1;
                double v2 = 2 * u2 - 1;

                double w = v1 * v1 + v2 * v2;

                while (w >= 1)
                {
                    u1 = random.NextDouble();
                    u2 = random.NextDouble();

                    v1 = 2 * u1 - 1;
                    v2 = 2 * u2 - 1;

                    w = v1 * v1 + v2 * v2;
                }

                double z1 = Math.Sqrt((-2 * Math.Log(w)) / w) * v1;
    
                values[i] = z1;
            }

            return values;

        }

        public static double[] RandomValues1(int count)
        {  // Borrowed from Project 4: generating a standard normal variable using polar rejection
            Random random = new Random();
            double[] values = new double[count];

            for (int i = 0; i < count; i++)
            {
                double u1 = random.NextDouble();
                double u2 = random.NextDouble();
                // transform uniform distribution interval to [-1,1]
                double v1 = 2 * u1 - 1;
                double v2 = 2 * u2 - 1;

                double w = v1 * v1 + v2 * v2;

                while (w >= 1)
                {
                    u1 = random.NextDouble();
                    u2 = random.NextDouble();

                    v1 = 2 * u1 - 1;
                    v2 = 2 * u2 - 1;

                    w = v1 * v1 + v2 * v2;
                }

                double z1 = Math.Sqrt((-2 * Math.Log(w)) / w) * v1;
                double z2 = Math.Sqrt((-2 * Math.Log(w)) / w) * v2;

                values[i] = z1;
            }

            return values;
    
        }
  
    }

}

  
