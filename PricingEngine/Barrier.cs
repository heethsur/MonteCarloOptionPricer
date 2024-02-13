using System;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Security.Cryptography.X509Certificates;

namespace MonteCarloOptionsPricer.PricingEngine
{
    public class Barrier : Option
    {
        public string? BarrierType {get ; set;}
        public double BarrierLevel {get; set;}
        public override double[] OptionPayoffs(double[,] paths)
        {// identify if option is valid for each path, and if so, apply European payoff
            double[] payoffs = new double[paths.GetLength(0)];
            double[] pathmin = PricePaths.GetPathMinimum(paths);
            double[] pathmax = PricePaths.GetPathMaximum(paths);
            if (BarrierType == "Down-Out") // if down and out
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
            else if (BarrierType == "Down-In") // if down and in
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
            else if (BarrierType == "Up-Out") // if up and out
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

}