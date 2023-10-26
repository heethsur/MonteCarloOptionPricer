using System;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Security.Cryptography.X509Certificates;

namespace MonteCarloOptionPricing
{
    class Program
    {
        static void Main(string[] args)
        {

            Stopwatch stopwatch = new(); // to measure implementation speed
           
            Console.WriteLine("Welcome to MonteCarlo Option Pricer!");
            Option option;

            // Ask user for parameters
            Console.WriteLine("Choose an option class (European, Asian, Digital, Barrier, Lookback)");
            string? optionClass = Console.ReadLine();
            switch (optionClass)
            {
                case "European":
                    option = new European();
                    break;
                case "Asian":
                    option = new Asian();
                    break;
                case "Digital":
                    option = new Digital();
                    Console.WriteLine("Enter payout amount: ");
                    ((Digital)option).PayoutAmount = Convert.ToDouble(Console.ReadLine());
                    break;
                case "Lookback":
                    option = new Lookback();
                    break;
                case "Barrier":
                    option = new Barrier();
                    Console.WriteLine("Enter Barrier Type: (DO for Down and Out, DI for Down and In, UO for Up and Out, UI for Up and In)");
                    ((Barrier)option).BarrierType = Console.ReadLine()?.ToUpper();
                    Console.WriteLine("Enter Barrier Level: ");
                    ((Barrier)option).BarrierLevel = Convert.ToDouble(Console.ReadLine());
                    break;
                default:
                    Console.WriteLine("Invalid input. ");
                    return;
            }   
            Console.WriteLine("Enter option type (true for Call, false for Put):");
            option.Call = Convert.ToBoolean(Console.ReadLine());           

            Console.WriteLine("Enter underlying price:");
            double Underlying = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter strike price:");
            option.Strike = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter time to expiration (in years):");
            option.TimeToExpiration = Convert.ToDouble(Console.ReadLine());

            Console.WriteLine("Enter volatility:");
            double Volatility = Convert.ToDouble(Console.ReadLine());
            option.Volatility = Volatility;

            Console.WriteLine("Enter risk-free rate:");
            double RiskFreeRate = Convert.ToDouble(Console.ReadLine());
            option.RiskFreeRate = RiskFreeRate;

            Console.WriteLine("Use Antithetic Sampling? (true/false)");
            bool antitheticsampling = Convert.ToBoolean(Console.ReadLine());

            Console.WriteLine("Use Control Variate? (true/false)");
            option.UseControlVariate = Convert.ToBoolean(Console.ReadLine());

            Console.WriteLine("Use MultiThreading? (true/false)");
            PricePaths.MultiThreading = Convert.ToBoolean(Console.ReadLine());     

            // instantiate normal random variables and paths. Apply payoffs, generate prices and greeks
            stopwatch.Start();
            Console.WriteLine("Generating random values...");
            double[,] normalrvs = RandomValues.NormalRandomValues(1000000, 12, antitheticsampling); // 1 million simulations with 12 timesteps
            stopwatch.Stop();
            double normtime = stopwatch.ElapsedMilliseconds;
            stopwatch.Reset();
            Console.WriteLine("Random values generated! (" + normtime/1000 + "s)");
            stopwatch.Start();
            Console.WriteLine("Generating price paths...");
            double[,] pricepaths = PricePaths.GeneratePricePaths(normalrvs, RiskFreeRate, Volatility, Underlying, option.TimeToExpiration);
            double[,] pricepaths1 = PricePaths.GeneratePricePaths(normalrvs, RiskFreeRate+0.01, Volatility, Underlying, option.TimeToExpiration);
            double[,] pricepaths2 = PricePaths.GeneratePricePaths(normalrvs, RiskFreeRate, Volatility+0.01, Underlying, option.TimeToExpiration);
            double[,] pricepaths3 = PricePaths.GeneratePricePaths(normalrvs, RiskFreeRate, Volatility, Underlying+1, option.TimeToExpiration);
            double[,] pricepaths4 = PricePaths.GeneratePricePaths(normalrvs, RiskFreeRate, Volatility, Underlying+2, option.TimeToExpiration);
            double[,] pricepaths5 = PricePaths.GeneratePricePaths(normalrvs, RiskFreeRate, Volatility, Underlying, option.TimeToExpiration + 0.01);
            stopwatch.Stop();
            double pathtime = stopwatch.ElapsedMilliseconds;
            stopwatch.Reset();
            Console.WriteLine("Price paths generated! (" + pathtime/1000 + "s)");
            stopwatch.Start();
            Console.WriteLine("Calculating payoffs...");
            double[] payoffs = option.Payoffs(pricepaths);
            double[] payoffs1 = option.Payoffs(pricepaths1);
            double[] payoffs2 = option.Payoffs(pricepaths2);
            double[] payoffs3 = option.Payoffs(pricepaths3);
            double[] payoffs4 = option.Payoffs(pricepaths4);
            double[] payoffs5 = option.Payoffs(pricepaths5);
            stopwatch.Stop();
            double payofftime = stopwatch.ElapsedMilliseconds;
            stopwatch.Reset();
            Console.WriteLine("Payoffs calculated! (" + payofftime/1000 + "s)");
            stopwatch.Start();
            Console.WriteLine("Calculating option price and greeks...");
            double Price = option.OptionPrice(payoffs, RiskFreeRate);
            double StdError = Option.StdError(payoffs);
            double Delta = option.OptionPrice(payoffs3, RiskFreeRate) - Price;
            double Gamma = option.OptionPrice(payoffs4, RiskFreeRate) - option.OptionPrice(payoffs3, RiskFreeRate) - Delta;
            double Vega = (option.OptionPrice(payoffs2, RiskFreeRate) - Price)*100;
            double Rho = (option.OptionPrice(payoffs1, RiskFreeRate+0.01) - Price)*100;
            option.TimeToExpiration += 0.01;
            double Theta = (option.OptionPrice(payoffs5, RiskFreeRate) - Price)*-100;
            stopwatch.Stop();
            double pricetime = stopwatch.ElapsedMilliseconds;
            stopwatch.Reset();
            Console.WriteLine("Prices and greeks calculated! (" + pricetime/1000 + "s)");
            
            // output results           
            Console.WriteLine("----------- "); 
            Console.WriteLine("Option price: " + Price);            
            Console.WriteLine("Standard error: " + StdError);            
            Console.WriteLine("Delta: " + Delta);           
            Console.WriteLine("Gamma: " + Gamma);     
            Console.WriteLine("Vega: " + Vega);           
            Console.WriteLine("Rho: " + Rho);         
            Console.WriteLine("Theta: " + Theta);
               


        }

        
    }

    public class Option
    {// a parent Option class
        public bool Call {get ; set;}
        public double Strike {get; set;}
        public double TimeToExpiration {get; set;}
        public double RiskFreeRate {get; set;}
        public double Volatility {get; set;}
        public bool UseControlVariate {get; set;}
        
        public virtual double[] OptionPayoffs(double[,] paths)
        {// default payoff generator
            double[] payoffs = new double[paths.GetLength(0)];
            return payoffs;
        }

        public virtual double[] CVPayoffs(double[,] paths)
        {// default control variate is underlying
            double[] payoffs = new double[paths.GetLength(0)];
            for (int i = 0; i < paths.GetLength(0); i++)
            {
                payoffs[i] = paths[i,paths.GetLength(1)-1];
            }
            return payoffs;
        }

        private static double CVBeta(double[] option_payoffs, double[] cv_payoffs)
        {// beta = cov(option_payoffs, controlvariate_payoffs)/var(controlvariate_payoffs)
            double covariance = 0;
            double variance_cv = 0;
            double beta = 0;
            double o_avg = option_payoffs.Average();
            double c_avg = cv_payoffs.Average();
            for(int i = 0; i < option_payoffs.Length; i++)
            {
                covariance += (option_payoffs[i] - o_avg)*(cv_payoffs[i] - c_avg);
                variance_cv += (cv_payoffs[i] - c_avg)*(cv_payoffs[i] - c_avg);
            }
            covariance /= option_payoffs.Length;
            variance_cv /= option_payoffs.Length;
            beta = covariance/variance_cv;
            return beta;
        }

        public virtual double[] Payoffs(double[,] paths)
        {// generate variance minimized payoffs(if controlvariate) otherwise return normal option payoffs
            double[] payoffs = new double[paths.GetLength(0)];
            double[] option_payoffs = OptionPayoffs(paths);
            double beta = 0;
            double[] cv_payoffs = new double[paths.GetLength(0)];
            if (UseControlVariate)
            {
                cv_payoffs = CVPayoffs(paths);
                beta = CVBeta(option_payoffs, cv_payoffs);              
                
            }
            double c_avg =  cv_payoffs.Average();
            for(int i = 0; i < option_payoffs.Length; i++)
            {
                payoffs[i] = option_payoffs[i] - (UseControlVariate? beta*(cv_payoffs[i] -  c_avg) : 0); 
            }
            return payoffs;
        }

       
        public double OptionPrice(double[] payoffs, double r)
        {// average and discount payoffs
            double discounted_payoff = payoffs.Average()*Math.Exp(-r*TimeToExpiration);
            return discounted_payoff;
        }

        public static double StdError(double[] payoffs)
        {// calculate std error of payoffs
            double[] sum_squared_payoffs = new double[payoffs.Length];
            for (int i = 0; i < payoffs.Length; i++)
            {
                sum_squared_payoffs[i] = payoffs[i] * payoffs[i];            
            }
            double avg_payoff =  payoffs.Average();
            double avg_squared_payoff = sum_squared_payoffs.Average();
            double stderror = Math.Sqrt(avg_squared_payoff - (avg_payoff * avg_payoff))/payoffs.Length;
            return stderror;
        }

    }

    public class Barrier : Option
    {
        public string? BarrierType {get ; set;}
        public double BarrierLevel {get; set;}
        public override double[] OptionPayoffs(double[,] paths)
        {// identify if option is valid for each path, and if so, apply European payoff
            double[] payoffs = new double[paths.GetLength(0)];
            double[] pathmin = PricePaths.GetPathMinimum(paths);
            double[] pathmax = PricePaths.GetPathMaximum(paths);
            if (BarrierType == "DO") // if down and out
            {
                for(int i = 0; i < paths.GetLength(0); i++) // for every path
                {
                    if (pathmin[i] > BarrierLevel) // path hasn't breached lower barrier?
                    {
                        if (Call) // is call?
                        {
                            payoffs[i] = Math.Max(paths[i, paths.GetLength(1) - 1] - Strike, 0);
                        }
                        else // is put?
                        {
                            payoffs[i] = Math.Max(Strike - paths[i, paths.GetLength(1) - 1], 0);
                        }
                    }
                    else // path has breached lower barrier?
                    {
                        payoffs[i] = 0;
                    }
                }
                return payoffs;
            }
            else if (BarrierType == "DI") // if down and in
            {
                for(int i = 0; i < paths.GetLength(0); i++) // for every path
                {
                    if (pathmin[i] < BarrierLevel) // path has breached lower barrier?
                    {
                        if (Call) // is call?
                        {
                            payoffs[i] = Math.Max(paths[i, paths.GetLength(1) - 1] - Strike, 0);
                        }
                        else // is put?
                        {
                            payoffs[i] = Math.Max(Strike - paths[i, paths.GetLength(1) - 1], 0);
                        }
                    }
                    else // path hasn't breached lower barrier?
                    {
                        payoffs[i] = 0;
                    }
                }
                return payoffs;
            }
            else if (BarrierType == "UO") // if up and out
            {
                for(int i = 0; i < paths.GetLength(0); i++) // for every path
                {
                    if (pathmax[i] < BarrierLevel) // path hasn't breached upper barrier?
                    {
                        if (Call) // is call?
                        {
                            payoffs[i] = Math.Max(paths[i, paths.GetLength(1) - 1] - Strike, 0);
                        }
                        else // is put?
                        {
                            payoffs[i] = Math.Max(Strike - paths[i, paths.GetLength(1) - 1], 0);
                        }
                    }
                    else // path has breached upper barrier?
                    {
                        payoffs[i] = 0;
                    }
                }
                return payoffs;
            }
            else // if up and in
            {
                for(int i = 0; i < paths.GetLength(0); i++) // for every path
                {
                    if (pathmax[i] > BarrierLevel) // path has breached upper barrier?
                    {
                        if (Call) // is call?
                        {
                            payoffs[i] = Math.Max(paths[i, paths.GetLength(1) - 1] - Strike, 0);
                        }
                        else // is put?
                        {
                            payoffs[i] = Math.Max(Strike - paths[i, paths.GetLength(1) - 1], 0);
                        }
                    }
                    else // path hasn't breached upper barrier?
                    {
                        payoffs[i] = 0;
                    }
                }
                return payoffs;
            }

        }
     

    }

    public class Lookback : Option
    {
        public override double[] OptionPayoffs(double[,] paths)
        {// get maximum price from each path and calculate payoff based on put or call
            double[] payoffs = new double[paths.GetLength(0)];
            double[] pathmax = PricePaths.GetPathMaximum(paths);
            for (int i = 0; i < paths.GetLength(0); i++)
            {
                if (Call) 
                { // if call
                    payoffs[i] = Math.Max(pathmax[i] - Strike, 0);
                }
                else
                { // else put
                    payoffs[i] = Math.Max(Strike - pathmax[i], 0);
                }
            }
            return payoffs;
        }  

    }

    public class Digital : Option
    {
        public double PayoutAmount {get; set;}
        public override double[] OptionPayoffs(double[,] paths)
        {// get terminal price from each path and calculate payoff based on payout amount and put or call
            double[] payoffs = new double[paths.GetLength(0)];
            for (int i = 0; i < paths.GetLength(0); i++)
            {
                if (Call) 
                { // if call
                    payoffs[i] = Math.Max(paths[i, paths.GetLength(1) - 1] - Strike, 0)/(paths[i, paths.GetLength(1) - 1] - Strike) * PayoutAmount;
                }
                else
                { // else put
                    payoffs[i] = Math.Max(Strike - paths[i, paths.GetLength(1) - 1], 0)/(Strike - paths[i, paths.GetLength(1) - 1]) * PayoutAmount;
                }
            }
            return payoffs;
        }
    }

    public class Asian : Option
    {
        public override double[] OptionPayoffs(double[,] paths)
        {// get average price from each path and calculate payoff based on put or call
            double[] payoffs = new double[paths.GetLength(0)];
            double[] pathavg = PricePaths.GetPathAverage(paths);
            for (int i = 0; i < paths.GetLength(0); i++)
            {
                if (Call) 
                { // if call
                    payoffs[i] = Math.Max(pathavg[i] - Strike, 0);
                }
                else
                { // else put
                    payoffs[i] = Math.Max(Strike - pathavg[i], 0);
                }
            }
            return payoffs;
        }   

    }

    public class European : Option
    {
        public override double[] OptionPayoffs(double[,] paths)
        {// get terminal price from each path and calculate payoff based on put or call
            double[] payoffs = new double[paths.GetLength(0)];
            for (int i = 0; i < paths.GetLength(0); i++)
            {
                if (Call) 
                { // if call
                    payoffs[i] = Math.Max(paths[i, paths.GetLength(1) - 1] - Strike, 0);
                }
                else
                { // else put
                    payoffs[i] = Math.Max(Strike - paths[i, paths.GetLength(1) - 1], 0);
                }
            }
            return payoffs;
        }

         public override double[] CVPayoffs(double[,] paths)
         {// Use delta-hedge as control variate
             double[] totalpnls = new double[paths.GetLength(0)]; //instantiate total pnls
             for(int i = 0; i < paths.GetLength(0); i++) // for each path
             {
                double totalpnl = 0; //initialize total pnl for each path
                for(int j = 1; j < paths.GetLength(1); j++) // for each timestep
                {
                    double t2x_prev = TimeToExpiration - (j-1)*TimeToExpiration/(paths.GetLength(1) - 1); // determine time to expiration as of previous period
                    double t2x_curr = TimeToExpiration - j*TimeToExpiration/(paths.GetLength(1) - 1); // determine time to expiration as of current period
                    double option_price_prev = AncillaryMethods.BSMPrice(paths[i,j-1],Strike,t2x_prev,RiskFreeRate,Volatility,Call); //option value at previous period
                    double option_price_curr = AncillaryMethods.BSMPrice(paths[i,j],Strike,t2x_curr,RiskFreeRate,Volatility,Call); //option value currently
                    double option_pnl = option_price_curr - option_price_prev; //pnl from holding the option between previous period to current period
                    double delta = AncillaryMethods.Delta(paths[i,j-1],Strike,t2x_prev,RiskFreeRate,Volatility,Call); // option delta at previous period
                    double underlying_pnl = -delta * (paths[i,j] - paths[i,j-1]); // hold negative delta amount of underlying, and calculate pnl from holding underlying
                    totalpnl += option_pnl + underlying_pnl; // calculate total pnl from holding underlying and option
                }
                totalpnls[i] = totalpnl;
             }
             return totalpnls;
         }

        
    }

    public class PricePaths
    {

        public static bool MultiThreading { get; set; } // using muti-threading?
        public static double[,] GeneratePricePaths(double[,] epsilon, double mu, double sigma, double S, double t2x)
        { // Generate price paths using Geometric Brownian Motion
            int numCores = Environment.ProcessorCount; // no. of cores
            int pathsperthread = epsilon.GetLength(0)/(numCores);
            double dt = t2x/epsilon.GetLength(1); // initialize timestep
            double[,] paths = new double[epsilon.GetLength(0), epsilon.GetLength(1) + 1]; // instantiate paths array
            if (MultiThreading) // if multithreading
            {
                Parallel.For(0, numCores, x => // for each thread
                {
                    int startindex = x * pathsperthread; //determine index of simulation that each thread begins execution at
                    int endindex = (x == (numCores - 1))? epsilon.GetLength(0): startindex + pathsperthread; //determine end index for each thread
                    for (int y = startindex; y < endindex; y++) //for each simulation in thread
                    {
                        paths[y,0] = S; // assign price at t_0 as S
                        for (int z = 1; z < (epsilon.GetLength(1) + 1); z++) // for every time step starting at t_1
                        {
                            paths[y,z] = paths[y,z-1] * Math.Exp((mu - (sigma * sigma/2)) * dt + (sigma *Math.Sqrt(dt)* epsilon[y,z-1]));
                        }
                    }
                });
            }
            else
            {
                for (int i = 0; i < epsilon.GetLength(0); i++) // for every simulation
                {
                    paths[i,0] = S; // assign price at t_0 as S
                    for (int j = 1; j < (epsilon.GetLength(1) + 1); j++) // for every time step starting at t_1
                    {// apply Geometric brownian motion to price at t - 1
                        paths[i,j] = paths[i,j-1] * Math.Exp((mu - (sigma * sigma/2)) * dt + (sigma *Math.Sqrt(dt)* epsilon[i,j-1]));
                    }
                }
            }
            return paths;
        }

        public static double[] GetPathAverage(double[,] paths)
        { // Get average prices along each path  
            double[] avgPrices = new double[paths.GetLength(0)];
            for (int i = 0; i < paths.GetLength(0); i++)
            {
                double[] path = new double[paths.GetLength(1)];
                for (int j = 0; j < paths.GetLength(1); j++)
                {
                    path[j] = paths[i,j];
                }
                avgPrices[i] = path.Average();
            }
            return avgPrices;   
        }

        public static double[] GetPathMaximum(double[,] paths)
        { // Get maximum prices along each path
            double[] maxPrices = new double[paths.GetLength(0)];
            for (int i = 0; i < paths.GetLength(0); i++)
            {
                double[] path = new double[paths.GetLength(1)];
                for (int j = 0; j < paths.GetLength(1); j++)
                {
                    path[j] = paths[i,j];
                }
                maxPrices[i] = path.Max();
            }
            return maxPrices;   
        }

        public static double[] GetPathMinimum(double[,] paths)
        { // Get minimum prices along each path
            double[] minPrices = new double[paths.GetLength(0)];
            for (int i = 0; i < paths.GetLength(0); i++)
            {
                double[] path = new double[paths.GetLength(1)];
                for (int j = 0; j < paths.GetLength(1); j++)
                {
                    path[j] = paths[i,j];
                }
                minPrices[i] = path.Min();
            }
            return minPrices;   
        }
    }

    public class RandomValues
    { // This class contains method to generate normal random values

        private static readonly Random random = new();

        public static double[,] NormalRandomValues(int simulations, int timesteps, bool UseAntitheticSampling)
        { // Generate an array of Random Normal Variables of length count using Polar Rejection
           double[] uniform_sample = GenerateVanderCorputSequence(2,1000); // Generate Van Der Corput Sequence of length 1000 with base 2
           double[,] randn = new double[simulations * (UseAntitheticSampling? 2 : 1), timesteps]; // instantiate output array of size count (double if antithetic)
           for (int j = 0; j < simulations; j++)
           {
                for (int i = 0; i < timesteps; i++) // loop through every index of output array
                {
                    double u1 = uniform_sample[random.Next(0, uniform_sample.Length)]; // randomly sample a value from the Van der Corput Sequence
                    double u2 = uniform_sample[random.Next(0, uniform_sample.Length)]; // randomly sample a value from the Van der Corput Sequence
                    double v1 = 2 * u1 - 1; // rescale u1 to be between -1 to 1
                    double v2 = 2 * u2 - 1; // rescale u1 to be between -1 to 1
                    double s = v1 * v1 + u2 * u2; 
                    if (s > 1) // until s < 1, resample u1 and u2
                    {
                        while (s > 1)
                        {
                            u1 = uniform_sample[random.Next(0, uniform_sample.Length)];
                            u2 = uniform_sample[random.Next(0, uniform_sample.Length)];
                            v1 = 2 * u1 - 1;
                            v2 = 2 * u2 - 1;
                            s = v1 * v1 + u2 * u2;
                        }
                    }                
                    double n = v1 * Math.Sqrt((-2 * Math.Log(s))/s); // normalize s
                    randn[j,i] = n; // assign normalized s to output array
                    if (UseAntitheticSampling)
                    {
                        randn[j + simulations, i] = -n; // assign antithetic complement 
                    }
                }
           }           
           return randn;
        }
        private static double[] GenerateVanderCorputSequence(int base_value, int n)
        {// Generate Van der Corput Sequence given a base value
            double[] sequence = new double[n];
            for (int i = 0; i < sequence.Length; i++)
            {
                int denom = base_value;
                double value = 0;
                int current = i + 1;
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
    }

    public class AncillaryMethods
    {// Other static methods for arithmetics and BSM calculations
        
        private static double D1(double S, double K, double t2x, double r, double v)
        {
            return (Math.Log(S/K) + (r+(v*v/2))*t2x)/(v*Math.Sqrt(t2x));
        }

        public static double Delta(double S, double K, double t2x, double r, double v, bool isCall)
        {
            return isCall? Cdf(D1(S,K,t2x,r,v)) : Cdf(-D1(S,K,t2x,r,v));
        } 

        public static double BSMPrice(double S, double K, double t2x, double r, double v, bool isCall)
        {
            return isCall? S*Cdf(D1(S,K,t2x,r,v)) - K*Math.Exp(-r*t2x)*Cdf(D1(S,K,t2x,r,v) - v*Math.Sqrt(t2x)) : -S*Cdf(D1(S,K,t2x,r,v)) + K*Math.Exp(-r*t2x)*Cdf(-D1(S,K,t2x,r,v) + v*Math.Sqrt(t2x));
        }

        private static double Pdf(double x)
        {
            return Math.Exp(-0.5 * x * x) / Math.Sqrt(2 * 3.141592653589793238462643383279502884197);
        }

        private static double Cdf(double x)
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


    }

}

    