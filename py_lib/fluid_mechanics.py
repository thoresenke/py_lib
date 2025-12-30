import math

G = 9.81

def compute_fluid_props(temperature_c: float, salinity_ppt: float, fluid_type: str = "water", pressure_pa: float = 101325.0):
    """
    Calculate density (kg/m3) and dynamic viscosity (Pa*s) for water or air.
    """
    if fluid_type.lower() == "air":
        # Air properties
        temp_k = temperature_c + 273.15
        
        # Density: Ideal Gas Law: rho = P / (R * T)
        # R_specific_air = 287.05 J/(kg*K)
        r_specific = 287.05
        density = pressure_pa / (r_specific * temp_k)
        
        # Viscosity: Sutherland's Law
        # mu = mu_ref * (T/T_ref)^(3/2) * (T_ref + S) / (T + S)
        mu_ref = 1.716e-5
        t_ref = 273.15
        s_sutherland = 110.4
        
        mu = mu_ref * (temp_k / t_ref)**1.5 * (t_ref + s_sutherland) / (temp_k + s_sutherland)
        
        return density, mu
    else:
        # Water properties
        temp = max(-5, min(50, temperature_c))
        sal = max(0, min(40, salinity_ppt))
        
        # Density (approx): base freshwater minus thermal expansion + salinity contribution
        density = 999.84 - 0.0631 * temp - 0.008523 * temp**2 + 0.0000665 * temp**3 + sal * 0.8
        
        # Dynamic viscosity (Pa·s) via Vogel equation with slight salinity bump
        mu_pure = 2.414e-5 * 10**(247.8 / (temp + 133.15))
        mu = mu_pure * (1 + 0.002 * (sal / 35))
        
        return density, mu

def calculate_reynolds(flow_rate: float, diameter_m: float, density: float, viscosity: float) -> float:
    if flow_rate == 0: return 0.0
    area = math.pi * (diameter_m / 2)**2
    velocity = flow_rate / area
    return (density * velocity * diameter_m) / viscosity

def calculate_pressure_loss(flow_rate: float, diameter_m: float, roughness_m: float, length_m: float, density: float, viscosity: float) -> float:
    """
    Calculate pressure loss for a straight pipe in Pa.
    """
    if flow_rate == 0: return 0.0
    
    area = math.pi * (diameter_m / 2)**2
    velocity = flow_rate / area
    
    re = (density * velocity * diameter_m) / viscosity
    
    if re < 2000:
        f = 64 / re
    else:
        # Swamee-Jain equation approximation for Darcy friction factor
        term = (roughness_m / (3.7 * diameter_m)) + (5.74 / (re**0.9))
        f = 0.25 / (math.log10(term))**2
        
    head_loss = f * (length_m / diameter_m) * (velocity**2 / (2 * G))
    pressure_loss_pa = density * G * head_loss
    
    return pressure_loss_pa

def calculate_component_loss(flow_rate: float, config: dict, rho: float, mu: float) -> float:
    """
    Calculate pressure loss for a pipe component (straight, bend, choke) in Pa.
    config expects: type, diameter (m), roughness (m), length (m), angle, radius (m), orificeDiameter (m), cd
    """
    if flow_rate == 0: return 0.0
    
    dia_m = config.get('diameter', 0)
    area = math.pi * (dia_m / 2)**2
    velocity = flow_rate / area
    
    head_loss = 0.0
    comp_type = config.get('type', 'straight')
    
    if comp_type == 'straight':
        rough_m = config.get('roughness', 0)
        re = (rho * velocity * dia_m) / mu
        if re < 2000:
            f = 64 / re
        else:
            term = (rough_m / (3.7 * dia_m)) + (5.74 / (re**0.9))
            f = 0.25 / (math.log10(term))**2
        head_loss = f * (config.get('length', 0) / dia_m) * (velocity**2 / (2 * G))
        
    elif comp_type == 'bend':
        r = max(config.get('radius', 0), config.get('diameter', 0) / 2)
        ratio = config.get('diameter', 0) / r
        theta = config.get('angle', 90)
        k = (0.131 + 0.1632 * (ratio**3.5)) * ((theta / 90)**0.5)
        head_loss = k * (velocity**2 / (2 * G))
        
    elif comp_type == 'choke':
        d2_m = config.get('orificeDiameter', 0)
        area2 = math.pi * (d2_m / 2)**2
        beta = d2_m / dia_m
        
        if d2_m >= dia_m or d2_m <= 0: return 0.0
        
        term1 = flow_rate / (config.get('cd', 0.6) * area2)
        term2 = rho * (1 - beta**4)
        dp = (term1**2) * term2 / 2
        return dp # Return directly in Pa
        
    pressure_loss_pa = rho * G * head_loss
    return pressure_loss_pa

def convert_pressure(pa: float, unit: str, rho: float) -> float:
    if unit == 'Bar': return pa / 1e5
    elif unit == 'Psi': return pa * 0.000145038
    elif unit == 'm': return pa / (rho * G)
    elif unit == 'mBar': return pa / 100
    else: return pa

def calculate_degassing_loss(num_nozzles: int, diameter_m: float, flow_rate: float, rho_air: float, cd: float) -> float:
    if num_nozzles <= 0 or diameter_m <= 0 or cd <= 0:
        return 0.0
        
    a_nozzle = math.pi * (diameter_m / 2)**2
    a_total = num_nozzles * a_nozzle
    
    if a_total == 0: return 0.0
    
    dp_pa = 0.5 * rho_air * (flow_rate / (cd * a_total))**2
    return dp_pa
