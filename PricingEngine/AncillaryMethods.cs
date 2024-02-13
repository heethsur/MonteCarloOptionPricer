using System;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Security.Cryptography.X509Certificates;

namespace MonteCarloOptionsPricer.PricingEngine
{
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