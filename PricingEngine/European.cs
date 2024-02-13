using System;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Security.Cryptography.X509Certificates;

namespace MonteCarloOptionsPricer.PricingEngine
{
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

}