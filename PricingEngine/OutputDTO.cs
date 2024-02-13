namespace MonteCarloOptionsPricer.PricingEngine
{
    public class OutputDTO
    {
        public double Price { get; set; }
        public string StandardError { get; set; }
        public double Delta { get; set; }
        public double Gamma { get; set; }
        public double Vega { get; set; }
        public double Rho { get; set; }
        public double Theta { get; set; }
        
    }
}