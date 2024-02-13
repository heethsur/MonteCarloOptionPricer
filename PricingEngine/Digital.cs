using System;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Security.Cryptography.X509Certificates;

namespace MonteCarloOptionsPricer.PricingEngine
{
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

}