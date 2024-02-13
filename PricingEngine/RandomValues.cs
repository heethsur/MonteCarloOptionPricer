using System;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;
using System.Security.Cryptography.X509Certificates;

namespace MonteCarloOptionsPricer.PricingEngine
{
    public class RandomValues
    { // This class contains method to generate normal random values

        public static bool UseAntitheticSampling { get; set; }
        private static readonly Random random = new();

        public static double[,] NormalRandomValues(int simulations, int timesteps)
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

}