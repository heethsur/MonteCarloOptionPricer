namespace MonteCarloOptionsPricer.PricingEngine
{
    public class OptionInputDTO
    {
        public double S { get; set; } 
        public double K { get; set; } 
        public double T { get; set; } 
        public bool Call {get; set; }
        public double Volatility { get; set; }
        public double RiskFreeRate { get; set; }
        public bool UseAntithetic { get; set; }
        public bool UseControlVariate { get; set; }
        public bool UseMultiThreading { get; set; }
        public double PayoutAmount { get; set; }
        public double BarrierLevel { get; set; }
        public string BarrierType { get; set; }
        
    }

}
