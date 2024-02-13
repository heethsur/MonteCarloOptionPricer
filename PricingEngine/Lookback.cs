using System;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Security.Cryptography.X509Certificates;

namespace MonteCarloOptionsPricer.PricingEngine
{
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

}