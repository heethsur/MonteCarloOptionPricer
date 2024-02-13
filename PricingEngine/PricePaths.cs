using System;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Security.Cryptography.X509Certificates;

namespace MonteCarloOptionsPricer.PricingEngine
{
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
                double maxInPath = paths[i, 0];
                for (int j = 1; j < paths.GetLength(1); j++)
                {
                    if (paths[i, j] > maxInPath)
                    {
                        maxInPath = paths[i, j];
                    }
                }
                maxPrices[i] = maxInPath;
            }
            return maxPrices;
            // double[] maxPrices = new double[paths.GetLength(0)];
            // for (int i = 0; i < paths.GetLength(0); i++)
            // {
            //     double[] path = new double[paths.GetLength(1)];
            //     for (int j = 0; j < paths.GetLength(1); j++)
            //     {
            //         path[j] = paths[i,j];
            //     }
            //     maxPrices[i] = path.Max();
            // }
            // return maxPrices;   
        }

        public static double[] GetPathMinimum(double[,] paths)
        { // Get minimum prices along each path
        double[] minPrices = new double[paths.GetLength(0)];
            for (int i = 0; i < paths.GetLength(0); i++)
            {
                double minInPath = paths[i, 0];
                for (int j = 1; j < paths.GetLength(1); j++)
                {
                    if (paths[i, j] < minInPath)
                    {
                        minInPath = paths[i, j];
                    }
                }
                minPrices[i] = minInPath;
            }
            return minPrices;
            // double[] minPrices = new double[paths.GetLength(0)];
            // for (int i = 0; i < paths.GetLength(0); i++)
            // {
            //     double[] path = new double[paths.GetLength(1)];
            //     for (int j = 0; j < paths.GetLength(1); j++)
            //     {
            //         path[j] = paths[i,j];
            //     }
            //     minPrices[i] = path.Min();
            // }
            // return minPrices;   
        }
    }

}