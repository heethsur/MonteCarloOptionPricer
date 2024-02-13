using System;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Security.Cryptography.X509Certificates;

namespace MonteCarloOptionsPricer.PricingEngine
{
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

}