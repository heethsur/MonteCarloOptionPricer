using System;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Security.Cryptography.X509Certificates;

namespace MonteCarloOptionsPricer.PricingEngine
{
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

}