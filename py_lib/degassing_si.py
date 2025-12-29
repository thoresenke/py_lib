import math
import sys
import os

# Add the parent directory to path so we can import py_lib
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from py_lib.units import u

def calculate_equilibrium_co2(temp_k: float, salinity_ppt: float, pressure_pa: float = 101325.0) -> float:
    """
    Calculates the equilibrium CO2 concentration (C*) in kg/m^3.
    For theory, see: file:///c:/Dev/fishfarm/py_lib/degassing_theory.md#1-equilibrium-co-concentration-c
    """
    
    # Convert SI inputs back to empirical units for the formula
    temp_c = u.DegtoC(temp_k)
    pressure_atm = pressure_pa / u.atm
    
    # Atmospheric CO2 partial pressure (approx 420 ppm)
    # pCO2_atm = 420e-6 * pressure_atm
    
    # Simplified engineering fit for C* (mg/L)
    c_star_mg_l = 0.6 * math.exp(-0.03 * (temp_c - 15)) * (1 - 0.005 * (salinity_ppt - 34))
    
    c_star_mg_l = max(0.1, c_star_mg_l)

    # Convert mg/L to kg/m^3 (SI)
    # 1 mg/L = 1 g/m^3 = 1e-3 kg/m^3
    return c_star_mg_l * (u.mg / u.L)

def calculate_ph(co2_conc_si: float, alkalinity_si: float, temp_k: float, salinity_ppt: float) -> float:
    """
    Estimates pH based on CO2 concentration and Alkalinity.
    For theory, see: file:///c:/Dev/fishfarm/py_lib/degassing_theory.md#2-carbonate-system-and-ph
    """
    # Convert SI inputs to empirical units
    co2_mg_l = co2_conc_si / (u.mg / u.L)
    alkalinity_meq_l = alkalinity_si / (u.mmol / u.L) # Assuming alkalinity_si is in mol/m^3 (eq/m^3)
    temp_c = u.DegtoC(temp_k)

    if co2_mg_l <= 0.0001:
        return 8.5 # Fallback for zero CO2
        
    # Convert units
    # MW CO2 = 44.01 g/mol
    co2_mol_kg = (co2_mg_l / 44010) 
    alk_mol_kg = alkalinity_meq_l / 1000.0

    # pK1 approximation (Millero et al., seawater scale)
    pK1 = 6.0 - 0.01 * (temp_c - 15)
    
    # pH = pK1 + log10(Alk / CO2)
    try:
        ph = pK1 + math.log10(alk_mol_kg / co2_mol_kg)
        return max(0.0, min(14.0, ph))
    except ValueError:
        return 7.0

def calculate_co2_out(
    co2_in_si: float,
    q_air_si: float,
    q_water_si: float,
    temp_k: float,
    salinity_ppt: float,
    efficiency_factor: float = 0.5
) -> dict:
    """
    Calculates the outlet CO2 concentration after degassing.
    For theory, see: file:///c:/Dev/fishfarm/py_lib/degassing_theory.md#3-mass-transfer-model-degassing
    
    Args:
        co2_in_si: Inlet CO2 concentration (kg/m^3)
        q_air_si: Airflow rate (m^3/s)
        q_water_si: Water flow rate (m^3/s)
        temp_k: Water temperature (Kelvin)
        salinity_ppt: Salinity (parts per thousand)
        efficiency_factor: Empirical efficiency factor (k).

    Returns:
        Dictionary containing:
        - 'co2_out': Outlet CO2 (kg/m^3)
        - 'removal_efficiency': Percentage of CO2 removed
        - 'c_star': Equilibrium CO2 (kg/m^3)
        - 'g_l_ratio': Gas to Liquid ratio
    """
    
    # 1. Calculate Equilibrium CO2 (C*)
    c_star_si = calculate_equilibrium_co2(temp_k, salinity_ppt)
    
    # 2. Calculate Gas/Liquid Ratio
    if q_water_si <= 0:
        return {
            "co2_out": co2_in_si, 
            "removal_efficiency": 0.0,
            "c_star": c_star_si,
            "g_l_ratio": 0.0
        }
        
    g_l_ratio = q_air_si / q_water_si
    
    # 3. Mass Transfer Model (Standard Aeration Equation)
    # Cout = C* + (Cin - C*) * e^(-k * G/L)
    exponent = -efficiency_factor * g_l_ratio
    transfer_efficiency = math.exp(exponent)
    
    co2_out_si = c_star_si + (co2_in_si - c_star_si) * transfer_efficiency
    
    # Ensure we don't go below C* (physically impossible via aeration alone)
    co2_out_si = max(co2_out_si, c_star_si)
    
    # Calculate removal percentage relative to available surplus
    removal_pct = 0.0
    if co2_in_si > c_star_si:
        removal_pct = (co2_in_si - co2_out_si) / (co2_in_si - c_star_si) * 100.0

    return {
        "co2_out": co2_out_si,
        "removal_efficiency": round(removal_pct, 2),
        "c_star": c_star_si,
        "g_l_ratio": round(g_l_ratio, 4)
    }

if __name__ == "__main__":
    # Test parameters
    temp       = u.CtoK(15.0)       # K
    salinity   = 34.0   # ppt
    alkalinity = 2.5 * (u.mmol / u.L)  # mol/m^3 (SI for mEq/L)
    co2_in     = 20.0 * (u.mg / u.L)     # kg/m^3
    water_flow = 500.0 * u.lpm # m^3/s
    efficiency = 0.5  # k factor

    print("Testing Degassing Physics (SI Units):")
    print(f"Temp: {temp} K, Salinity: {salinity} ppt")
    print(f"Water Flow: {water_flow} m^3/s, CO2 In: {co2_in} kg/m^3")
    print("-" * 75)
    print(f"{'Airflow (L/min)':<20} | {'CO2 Out (mg/L)':<15} | {'pH Out':<10} | {'Efficiency %':<10}")
    print("-" * 75)

    # Simulate Airflow from 0 to 1000 L/min
    for airflow_lpm in range(0, 1001, 100):
        airflow_si = airflow_lpm * u.lpm
        
        result = calculate_co2_out(
            co2_in_si=co2_in,
            q_air_si=airflow_si,
            q_water_si=water_flow,
            temp_k=temp,
            salinity_ppt=salinity,
            efficiency_factor=efficiency
        )
        
        co2_out_si = result['co2_out']
        ph_out = calculate_ph(co2_out_si, alkalinity, temp, salinity)
        
        # Convert back for display
        co2_out_display = co2_out_si / (u.mg / u.L)
        
        print(f"{airflow_lpm:<20} | {co2_out_display:<15.2f} | {ph_out:<10.2f} | {result['removal_efficiency']:<10.1f}")
