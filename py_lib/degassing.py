import math

def calculate_equilibrium_co2(temp_c: float, salinity_ppt: float, pressure_atm: float = 1.0) -> float:
    """
    Calculates the equilibrium CO2 concentration (C*) in mg/L based on 
    Temperature (C) and Salinity (ppt).
    
    Uses a simplified Weiss (1974) approximation for solubility.
    """
    # Kelvin
    T_K = temp_c + 273.15
    
    # Henry's constant K0 (mol/kg/atm) approximation
    # This is a simplified engineering fit. 
    # Base solubility at 15C, 35ppt is approx 0.6-0.7 mg/L for 400ppm air.
    
    # Atmospheric CO2 partial pressure (approx 420 ppm)
    pCO2_atm = 420e-6 * pressure_atm
    
    # Solubility decreases with Temp and Salinity
    # Base reference: ~0.6 mg/L at 15C, 34ppt
    # Temp effect: approx -3% per degree C
    # Salinity effect: approx -0.5% per ppt
    
    c_star = 0.6 * math.exp(-0.03 * (temp_c - 15)) * (1 - 0.005 * (salinity_ppt - 34))
    
    return max(0.1, c_star)

def calculate_ph(co2_mg_l: float, alkalinity_meq_l: float, temp_c: float, salinity_ppt: float) -> float:
    """
    Estimates pH based on CO2 concentration and Alkalinity using 
    carbonate system equilibrium approximations.
    """
    if co2_mg_l <= 0.0001:
        return 8.5 # Fallback for zero CO2
        
    # Convert units
    # MW CO2 = 44.01 g/mol
    co2_mol_kg = (co2_mg_l / 44010) 
    alk_mol_kg = alkalinity_meq_l / 1000.0

    # pK1 approximation (Millero et al., seawater scale)
    # Simplified fit: pK1 ~ 5.9 at 15C, 35ppt
    pK1 = 6.0 - 0.01 * (temp_c - 15)
    
    # pH = pK1 - log10(CO2 / HCO3)
    # Assuming Alk ~ [HCO3-] (valid for pH < 8.5 approx)
    # pH = pK1 + log10(Alk / CO2)
    try:
        ph = pK1 + math.log10(alk_mol_kg / co2_mol_kg)
        return max(0.0, min(14.0, ph))
    except ValueError:
        return 7.0

def calculate_co2_out(
    co2_in_mg_l: float,
    q_air_l_min: float,
    q_water_l_min: float,
    temp_c: float,
    salinity_ppt: float,
    efficiency_factor: float = 0.5
) -> dict:
    """
    Calculates the outlet CO2 concentration and resulting pH after degassing.
    
    Args:
        co2_in_mg_l: Inlet CO2 concentration (mg/L)
        q_air_l_min: Airflow rate (L/min)
        q_water_l_min: Water flow rate (L/min)
        temp_c: Water temperature (Celsius)
        salinity_ppt: Salinity (parts per thousand)
        efficiency_factor: Empirical efficiency factor (k) for the specific aerator type.
                           Typical range 0.3 - 0.8 depending on diffuser type.

    Returns:
        Dictionary containing:
        - 'co2_out': Outlet CO2 (mg/L)
        - 'removal_efficiency': Percentage of CO2 removed
        - 'c_star': Equilibrium CO2 (mg/L)
        - 'g_l_ratio': Gas to Liquid ratio
    """
    
    # 1. Calculate Equilibrium CO2 (C*)
    c_star = calculate_equilibrium_co2(temp_c, salinity_ppt)
    
    # 2. Calculate Gas/Liquid Ratio
    if q_water_l_min <= 0:
        return {
            "co2_out": co2_in_mg_l, 
            "removal_efficiency": 0.0,
            "c_star": c_star,
            "g_l_ratio": 0.0
        }
        
    g_l_ratio = q_air_l_min / q_water_l_min
    
    # 3. Mass Transfer Model (Standard Aeration Equation)
    # Cout = C* + (Cin - C*) * e^(-k * G/L)
    exponent = -efficiency_factor * g_l_ratio
    transfer_efficiency = math.exp(exponent)
    
    co2_out = c_star + (co2_in_mg_l - c_star) * transfer_efficiency
    
    # Ensure we don't go below C* (physically impossible via aeration alone)
    co2_out = max(co2_out, c_star)
    
    # Calculate removal percentage relative to available surplus
    removal_pct = 0.0
    if co2_in_mg_l > c_star:
        removal_pct = (co2_in_mg_l - co2_out) / (co2_in_mg_l - c_star) * 100.0

    return {
        "co2_out": round(co2_out, 4),
        "removal_efficiency": round(removal_pct, 2),
        "c_star": round(c_star, 4),
        "g_l_ratio": round(g_l_ratio, 4)
    }

if __name__ == "__main__":
    # Test parameters
    temp = 15.0       # C
    salinity = 34.0   # ppt
    alkalinity = 2.5  # mEq/L
    co2_in = 20.0     # mg/L
    water_flow = 500.0 # L/min
    efficiency = 0.5  # k factor

    print("Testing Degassing Physics (Field Units):")
    print(f"Temp: {temp} C, Salinity: {salinity} ppt")
    print(f"Water Flow: {water_flow} L/min, CO2 In: {co2_in} mg/L")
    print("-" * 75)
    print(f"{'Airflow (L/min)':<20} | {'CO2 Out (mg/L)':<15} | {'pH Out':<10} | {'Efficiency %':<10}")
    print("-" * 75)

    # Simulate Airflow from 0 to 1000 L/min
    for airflow in range(0, 1001, 100):
        result = calculate_co2_out(
            co2_in_mg_l=co2_in,
            q_air_l_min=airflow,
            q_water_l_min=water_flow,
            temp_c=temp,
            salinity_ppt=salinity,
            efficiency_factor=efficiency
        )
        
        co2_out = result['co2_out']
        ph_out = calculate_ph(co2_out, alkalinity, temp, salinity)
        
        print(f"{airflow:<20} | {co2_out:<15.2f} | {ph_out:<10.2f} | {result['removal_efficiency']:<10.1f}")
