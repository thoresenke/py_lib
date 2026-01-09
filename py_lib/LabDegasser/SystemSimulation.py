import numpy as np
from dataclasses import dataclass, field, asdict
import json
from typing import Optional, Dict, Any
from scipy.optimize import root_scalar

# Try imports to handle both package and script usage
try:
    from .Chiller import ChillerParams, calculate_capacity_ua_implicit
    from ..matLib import u
except ImportError:
    import sys
    import os
    # Add current directory to path for Chiller
    sys.path.append(os.path.dirname(__file__))
    # Add parent directory to path for matLib
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    
    from Chiller import ChillerParams, calculate_capacity_ua_implicit
    from matLib import u

@dataclass
class AmbientParams:
    temperature_C: float = 20.0

@dataclass
class BlowerParams:
    flow_rate_m3hr: float = 300.0
    pressure_rise_mbar: float = 150.0 # Typical for aeration
    efficiency: float = 0.6
    
    @property
    def estimated_heat_W(self) -> float:
        # Power = Q * dP / eta
        # Q in m3/s, dP in Pa
        q_m3s = self.flow_rate_m3hr / 3600.0
        dp_pa = self.pressure_rise_mbar * 100.0
        power_W = (q_m3s * dp_pa) / self.efficiency
        return power_W

@dataclass
class DiffuserParams:
    num_holes: int = 100
    hole_diameter_mm: float = 2.0
    water_depth_m: float = 1.5
    discharge_coeff: float = 0.62

    def calculate_pressure_drop_pa(self, flow_rate_m3hr: float) -> float:
        """
        Calculates pressure drop across diffuser holes + hydrostatic head.
        """
        # Hydrostatic
        rho_water = 1000.0
        g = 9.81
        dp_hydro = rho_water * g * self.water_depth_m
        
        # Orifice drop
        # v = Q / A
        # dp = 0.5 * rho_air * v^2 / Cd^2
        # Note: Using air density at depth? Or average?
        # Simplified: Incompressible air assumption (inaccurate but sufficient for estimation)
        rho_air = 1.2 # kg/m3
        
        q_m3s = flow_rate_m3hr / 3600.0
        area_hole = np.pi * (self.hole_diameter_mm / 2000.0)**2
        total_area = self.num_holes * area_hole
        
        if total_area <= 0:
            return dp_hydro
            
        v = q_m3s / total_area
        
        # Check for choked flow (approximate limit)
        if v > 300:
            print(f"Warning: Calculated diffuser velocity {v:.1f} m/s is very high.")
            
        dp_orifice = 0.5 * rho_air * (v**2) / (self.discharge_coeff**2)
        
        return dp_hydro + dp_orifice

@dataclass
class TankParams:
    volume_m3: float = 2.0
    count: int = 2
    surface_area_per_tank_m2: float = 6.0 # Approx for 2m3 tank
    h_conv: float = 5.0 # W/m2K natural convection

@dataclass
class PipeParams:
    length_m: float = 5.0
    diameter_mm: float = 10.0
    emissivity: float = 0.85
    h_conv: float = 10.0 # W/m2K

@dataclass
class PumpParams:
    flow_rate_lpm: float = 10.0
    power_W: float = 100.0 # Heat input to fluid

@dataclass
class ChillerSystemParams:
    chiller_params: ChillerParams
    pel_W: float = 1250.0
    enabled: bool = True

@dataclass
class SystemConfig:
    ambient: AmbientParams
    blower: BlowerParams
    diffuser: DiffuserParams
    tanks: TankParams
    pipes: PipeParams
    pump: PumpParams
    chiller: ChillerSystemParams

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=4)
    
    @classmethod
    def from_json(cls, json_str: str) -> 'SystemConfig':
        data = json.loads(json_str)
        # Reconstruct nested objects
        data['ambient'] = AmbientParams(**data['ambient'])
        data['blower'] = BlowerParams(**data['blower'])
        data['diffuser'] = DiffuserParams(**data['diffuser'])
        data['tanks'] = TankParams(**data['tanks'])
        data['pipes'] = PipeParams(**data['pipes'])
        data['pump'] = PumpParams(**data['pump'])
        
        chiller_data = data['chiller']
        chiller_params_data = chiller_data.pop('chiller_params')
        chiller_data['chiller_params'] = ChillerParams(**chiller_params_data)
        data['chiller'] = ChillerSystemParams(**chiller_data)
        
        return cls(**data)

def estimate_steady_state_temperature(config: SystemConfig) -> float:
    """
    Calculates the steady state water temperature of the system.
    Returns Temperature in Celsius.
    """
    
    # Constants
    SIGMA = 5.67e-8 # Stefan-Boltzmann constant
    
    # Unpack parameters for easier access
    Ta_C = config.ambient.temperature_C
    Ta_K = u.C2K(Ta_C)
    
    # Heat Sources (assumed constant)
    Q_pump = config.pump.power_W
    
    # Calculate Blower Heat based on Diffuser restriction
    # If pressure is not fixed in blower params (or we want to override), calculate it
    # For this simulation, we assume the blower must overcome the diffuser pressure
    dp_diffuser_pa = config.diffuser.calculate_pressure_drop_pa(config.blower.flow_rate_m3hr)
    
    # Power = Q * dP / eta
    q_m3s = config.blower.flow_rate_m3hr / 3600.0
    blower_power_W = (q_m3s * dp_diffuser_pa) / config.blower.efficiency
    
    Q_blower = blower_power_W
    Q_in = Q_pump + Q_blower
    
    # print(f"Debug: Blower dP={dp_diffuser_pa:.0f} Pa, Power={blower_power_W:.0f} W")
    
    # Geometry for losses
    pipe_area = np.pi * (config.pipes.diameter_mm / 1000.0) * config.pipes.length_m
    tank_area = config.tanks.count * config.tanks.surface_area_per_tank_m2
    
    def net_heat_balance(Tw_C):
        Tw_K = u.C2K(Tw_C)
        
        # 1. Chiller Cooling (Method 2: UA Implicit)
        Q_chiller = 0.0
        if config.chiller.enabled:
            # calculate_capacity_ua_implicit(Tw, Ta, Pel, eta_carnot, UA_evap, UA_cond, Qev_max)
            # Note: Chiller function expects Kelvin
            Q_chiller = calculate_capacity_ua_implicit(
                Tw_K, 
                Ta_K, 
                config.chiller.pel_W, 
                config.chiller.chiller_params.eta_carnot,
                config.chiller.chiller_params.UA_evap,
                config.chiller.chiller_params.UA_cond,
                config.chiller.chiller_params.Qev_max
            )
            if np.isnan(Q_chiller):
                Q_chiller = 0.0 # Chiller failed or out of range
        
        # 2. Pipe Losses (Radiation + Convection)
        # Convection
        Q_pipe_conv = config.pipes.h_conv * pipe_area * (Tw_K - Ta_K)
        # Radiation
        Q_pipe_rad = config.pipes.emissivity * SIGMA * pipe_area * (Tw_K**4 - Ta_K**4)
        Q_pipe_loss = Q_pipe_conv + Q_pipe_rad
        
        # 3. Tank Losses (Convection only for simplicity, or add rad if needed)
        Q_tank_loss = config.tanks.h_conv * tank_area * (Tw_K - Ta_K)
        
        # Total Balance
        # Residual = Heat In - Heat Out
        # We want Residual = 0
        residual = Q_in - (Q_chiller + Q_pipe_loss + Q_tank_loss)
        # print(f"Tw={Tw_C:.2f}, Qin={Q_in:.1f}, Qchiller={Q_chiller:.1f}, Qloss={Q_pipe_loss+Q_tank_loss:.1f}, Res={residual:.1f}")
        return residual

    # Solve for Tw_C
    # Search range: -10 to 200 (Theoretical)
    try:
        # Check endpoints first to avoid blind brentq failure
        low, high = -10.0, 200.0
        r_low = net_heat_balance(low)
        r_high = net_heat_balance(high)
        
        if np.isnan(r_low) or np.isnan(r_high):
             print(f"Warning: Solver endpoints returned NaN. Low({low})={r_low}, High({high})={r_high}")
             return np.nan
             
        if r_low * r_high > 0:
             print(f"Warning: Solver endpoints have same sign. Low({low})={r_low}, High({high})={r_high}")
             return np.nan

        sol = root_scalar(net_heat_balance, bracket=[low, high], method='brentq')
        if sol.converged:
            return sol.root
        else:
            raise RuntimeError("Solver did not converge")
    except ValueError:
        # Bracket error, try finding a valid bracket or just return nan
        return np.nan

if __name__ == "__main__":
    # Load configuration from JSON file
    config_path = os.path.join(os.path.dirname(__file__), 'system_config.json')
    
    try:
        with open(config_path, 'r') as f:
            json_str = f.read()
            sys_config = SystemConfig.from_json(json_str)
            print(f"Loaded configuration from {config_path}")
    except FileNotFoundError:
        print(f"Configuration file not found at {config_path}. Using defaults.")
        # Create default chiller params (similar to Titan 4000 from Chiller.py)
        chiller_p = ChillerParams(
            UA_evap=1500.0, # Placeholder
            UA_cond=3000.0, # Placeholder
            eta_carnot=0.4,
            lift_min=10.0,
            Qev_max=5000.0
        )
        
        sys_config = SystemConfig(
            ambient=AmbientParams(temperature_C=20.0),
            blower=BlowerParams(flow_rate_m3hr=300.0),
            diffuser=DiffuserParams(hole_diameter_mm=5.0),
            tanks=TankParams(volume_m3=2.0, count=2),
            pipes=PipeParams(length_m=5.0, diameter_mm=10.0),
            pump=PumpParams(flow_rate_lpm=10.0, power_W=100.0),
            chiller=ChillerSystemParams(chiller_params=chiller_p, pel_W=1250.0)
        )
    
    # print("System Configuration:")
    # print(sys_config.to_json())
    
    
    # Calculate and print blower stats for verification
    dp = sys_config.diffuser.calculate_pressure_drop_pa(sys_config.blower.flow_rate_m3hr)
    q_m3s = sys_config.blower.flow_rate_m3hr / 3600.0
    power = (q_m3s * dp) / sys_config.blower.efficiency
    print(f"\nCalculated Blower dP: {dp/100:.1f} mbar")
    print(f"Calculated Blower Heat Input: {power:.1f} W")
    
    T_steady = estimate_steady_state_temperature(sys_config)
    print(f"\nEstimated Steady State Temperature: {T_steady:.2f} °C")
