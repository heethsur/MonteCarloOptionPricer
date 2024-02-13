using Microsoft.AspNetCore.Mvc;

namespace MonteCarloOptionsPricer.Controllers;

[ApiController]
[Route("[controller]")]
public class OptionPriceController : ControllerBase
{
   [HttpPost("{option_class}")] // specify option class in URL. For example: OptionPrice/European, OptionPrice/Asian, OptionPrice/Barrier
   public ActionResult<PricingEngine.OutputDTO> CalculatePrice([FromBody] PricingEngine.OptionInputDTO input, string option_class)
   {// get option parameters from body JSON. Refer to OptionInputDTO for properties needed.
    // Example body input: {"S": 100,"K": 100,"T": 2,"Call": true,"Volatility": 0.25,"RiskFreeRate": 0.07,"UseAntithetic": true,"UseControlVariate": false,"UseMultiThreading": true,"PayoutAmount": 0,"BarrierLevel": 0, "BarrierType": "None"}
    // Note that values need to be entered into PayoutAmount, BarrierType, and BarrierLevel irrespective of option class. You can add dummy values if option class isn't Barrier or Digital.
        
        PricingEngine.OutputDTO result = new PricingEngine.OutputDTO(); // instantiate return class DTO
        PricingEngine.Option option = null; // instantiate option
        switch (option_class) // determine option class and instantiate option accordingly
            {
                case "European":
                    option = new PricingEngine.European();
                    break;
                case "Asian":
                    option = new PricingEngine.Asian();
                    break;
                case "Digital":
                    option = new PricingEngine.Digital();
                    ((PricingEngine.Digital)option).PayoutAmount = input.PayoutAmount;
                    break;
                case "Lookback":
                    option = new PricingEngine.Lookback();
                    break;
                case "Barrier":
                    option = new PricingEngine.Barrier();
                    ((PricingEngine.Barrier)option).BarrierType = input.BarrierType;
                    ((PricingEngine.Barrier)option).BarrierLevel = input.BarrierLevel;
                    break;
                default:
                    return BadRequest("Invalid option_class");
            }
            
        // assign option properties
        option.Call = input.Call;
        option.Strike = input.K;
        option.TimeToExpiration = input.T;
        option.RiskFreeRate = input.RiskFreeRate;
        option.Volatility = input.Volatility;
        option.UseControlVariate = input.UseControlVariate;
        // assign other properties       
        PricingEngine.RandomValues.UseAntitheticSampling = input.UseAntithetic;
        PricingEngine.PricePaths.MultiThreading = input.UseMultiThreading;
        // generate random values
        double[,] normalrvs = PricingEngine.RandomValues.NormalRandomValues(1000000, 12);
        // generate price paths for price and greeks
        double[,] pricepaths = PricingEngine.PricePaths.GeneratePricePaths(normalrvs, input.RiskFreeRate, input.Volatility, input.K, input.T);
        double[,] pricepaths1 = PricingEngine.PricePaths.GeneratePricePaths(normalrvs, input.RiskFreeRate+0.01, input.Volatility, input.K, input.T);
        double[,] pricepaths2 = PricingEngine.PricePaths.GeneratePricePaths(normalrvs, input.RiskFreeRate, input.Volatility+0.01, input.K, input.T);
        double[,] pricepaths3 = PricingEngine.PricePaths.GeneratePricePaths(normalrvs, input.RiskFreeRate, input.Volatility, input.K+1, input.T);
        double[,] pricepaths4 = PricingEngine.PricePaths.GeneratePricePaths(normalrvs, input.RiskFreeRate, input.Volatility, input.K+2, input.T);
        double[,] pricepaths5 = PricingEngine.PricePaths.GeneratePricePaths(normalrvs, input.RiskFreeRate, input.Volatility, input.K, input.T + 0.01);
        // generate payoffs for each price path
        double[] payoffs = option.Payoffs(pricepaths);
        double[] payoffs1 = option.Payoffs(pricepaths1);
        double[] payoffs2 = option.Payoffs(pricepaths2);
        double[] payoffs3 = option.Payoffs(pricepaths3);
        double[] payoffs4 = option.Payoffs(pricepaths4);
        double[] payoffs5 = option.Payoffs(pricepaths5);
        // calculate price, std error, and greeks. Add to output class DTO
        result.Price = Math.Round(option.OptionPrice(payoffs, option.RiskFreeRate),3);
        result.StandardError = PricingEngine.Option.StdError(payoffs).ToString("E3");
        result.Delta = Math.Round(option.OptionPrice(payoffs3, option.RiskFreeRate) - result.Price,3);
        result.Gamma = Math.Round(option.OptionPrice(payoffs4, option.RiskFreeRate) - option.OptionPrice(payoffs3, option.RiskFreeRate) - result.Delta,3);
        result.Vega = Math.Round((option.OptionPrice(payoffs2, option.RiskFreeRate) - result.Price)*100,3);
        result.Rho = Math.Round((option.OptionPrice(payoffs1, option.RiskFreeRate+0.01) - result.Price)*100,3);
        option.TimeToExpiration += 0.01;
        result.Theta = Math.Round((option.OptionPrice(payoffs5, option.RiskFreeRate) - result.Price)*-100,3);
        
        return Ok(result);
   }

    
}
