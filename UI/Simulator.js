

// dynamically hide or show additional option parameters if option type is Barrier or Digital

document.getElementById("OptionType").addEventListener("change", function() {
    if(this.value === "Digital") { // if digital, show box to input payout amount
        document.getElementById("Plabel").style.display = "block";
        document.getElementById("PayoutAmount").style.display = "block";
        document.getElementById("BTypelabel").style.display = "none";
        document.getElementById("BarrierType").style.display = "none";
        document.getElementById("BLevellabel").style.display = "none";
        document.getElementById("BarrierLevel").style.display = "none";
    } else if(this.value === "Barrier") { // if Barrier, show boxes to input barrier level and barrier type
        document.getElementById("BTypelabel").style.display = "block";
        document.getElementById("BarrierType").style.display = "block";
        document.getElementById("BLevellabel").style.display = "block";
        document.getElementById("BarrierLevel").style.display = "block";
        document.getElementById("Plabel").style.display = "none";
        document.getElementById("PayoutAmount").style.display = "none";
    } else { // otherwise keep them hidden
        document.getElementById("Plabel").style.display = "none";
        document.getElementById("PayoutAmount").style.display = "none";
        document.getElementById("BTypelabel").style.display = "none";
        document.getElementById("BarrierType").style.display = "none";
        document.getElementById("BLevellabel").style.display = "none";
        document.getElementById("BarrierLevel").style.display = "none";
    }
});

// make call to pricing API when clicking simulate button

document.getElementById("simulate").addEventListener("click", function(event) {

    let optionClass = document.getElementById("OptionType").value; // determine option class for URL
    let expdate = new Date(document.getElementById("End").value); 
    let startdate = new Date(document.getElementById("Start").value);
    let t2x = (expdate - startdate)/(1000*60*60*24*365); // calculate time to expiration
    let pa = "100"; // set default values for payout amount, barrier type and level since our input DTO needs them 
    let bt = "Down-Out"
    let bl = "100";
    if (optionClass === "Barrier") { // change barrier level and type from default to user input if barrier option is selected
        bt = document.getElementById("BarrierType").value;
        bl = document.getElementById("BarrierLevel").value;
    }
    if (optionClass === "Digital") { // change payout amount from default to user input if digi is selected
        pa = document.getElementById("PayoutAmount").value;
    }
    let inputData = { // populate inputDTO
        S: document.getElementById("S").value,
        K: document.getElementById("K").value,
        T: t2x,
        Call: document.getElementById("call").checked,
        Volatility: document.getElementById("V").value,
        RiskFreeRate: document.getElementById("R").value,
        UseAntiThetic: document.getElementById("Antithetic").checked,
        UseControlVariate: document.getElementById("ControlVariate").checked,
        UseMultiThreading: document.getElementById("MultiThread").checked,
        PayoutAmount: pa,
        BarrierLevel: bl,
        BarrierType: bt
    };

    console.log(JSON.stringify(inputData)); // to check if input looks okay

    let xhr = new XMLHttpRequest();
    xhr.open("POST",`http://localhost:5225/OptionPrice/${optionClass}`, true); // makw POST API request
    xhr.setRequestHeader("Content-Type", "application/json");
    xhr.onreadystatechange = () => {
        if(xhr.readyState === XMLHttpRequest.DONE){
            if(xhr.status === 200){
               console.log("Processing...")
               let result = JSON.parse(xhr.responseText);
               document.getElementById("Price").textContent = result.price; // extract results from outputDTO and populate table in UI
               document.getElementById("Delta").textContent = result.delta;
               document.getElementById("Gamma").textContent = result.gamma;
               document.getElementById("Vega").textContent = result.vega;
               document.getElementById("Rho").textContent = result.rho;
               document.getElementById("Theta").textContent = result.theta;
               document.getElementById("StdError").textContent = result.standardError;
               console.log("Done.")
            } else {
                console.error('Error in request:',xhr.responseText);
            }

        }
    };

    xhr.send(JSON.stringify(inputData)); 

});

// reset input and output values when clicking refresh button
document.getElementById("refresh").addEventListener("click", function(event) {
    document.getElementById("OptionType").value = "European"; 
    document.getElementById("End").value = "";
    document.getElementById("Start").value = "";
    document.getElementById("S").value = "";
    document.getElementById("K").value = "";
    document.getElementById("V").value = "";
    document.getElementById("R").value = "";
    document.getElementById("call").checked = false; 
    document.getElementById("put").checked = false;
    document.getElementById("Antithetic").checked = false;
    document.getElementById("ControlVariate").checked = false;
    document.getElementById("MultiThread").checked = false;
    document.getElementById("Price").textContent = "?";
    document.getElementById("Delta").textContent = "?";
    document.getElementById("Gamma").textContent = "?";
    document.getElementById("Vega").textContent = "?";
    document.getElementById("Rho").textContent = "?";
    document.getElementById("Theta").textContent = "?";
    document.getElementById("StdError").textContent = "?";
});
