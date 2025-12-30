def calculate_heat_input(fluid_type: str, flow_rate: float, temp_source: float, temp_target: float) -> float:
    """
    Calculate heat input in Watts.
    """
    fluids = {
        'air': {'rho': 1.225, 'cp': 1005},
        'water': {'rho': 1000, 'cp': 4186}
    }
    
    props = fluids.get(fluid_type, fluids['air'])
    rho = props['rho']
    cp = props['cp']
    
    delta_t = temp_source - temp_target
    power_watts = rho * flow_rate * cp * delta_t
    return power_watts
