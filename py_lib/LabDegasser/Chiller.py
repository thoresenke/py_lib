import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from scipy.optimize import root, brentq
try:
    from matLib import mCol, u
except ImportError:
    # Fallback if running as script from different dir without path setup
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from matLib import mCol, u

CP_WATER = 4186.0  # J/(kg*K)

@dataclass
class ChillerParams:
    """
    Parameters defining the physical characteristics and limits of the chiller.
    
    UA_evap: Evaporator heat transfer capability (Watts/Kelvin).(hA equivalent)
    UA_cond: Condenser heat transfer capability (Watts/Kelvin).
    eta_carnot:  Efficiency factor relative to ideal Carnot cycle (0.0 - 1.0).
    lift_min:  Minimum temperature lift (Condenser T - Evaporator T) to maintain operation.
    Tevap_min: Minimum allowed evaporator temperature (freeze protection).
    Tcond_max: Maximum allowed condenser temperature (safety limit).
    Qev_max:   Maximum cooling capacity (soft limit).
    """
    UA_evap: float
    UA_cond: float
    eta_carnot: float
    lift_min: float = 10.0
    Tevap_min: float = u.C2K(0.5)
    Tcond_max: float = u.C2K(70.0)
    Qev_max: float = 5000.0
    
    def __post_init__(self):
        if self.UA_evap <= 0:
            raise ValueError(f"UA_evap must be positive, got {self.UA_evap}")
        if self.UA_cond <= 0:
            raise ValueError(f"UA_cond must be positive, got {self.UA_cond}")
        if not 0 < self.eta_carnot <= 1:
            raise ValueError(f"eta_carnot must be in (0, 1], got {self.eta_carnot}")

def cop_carnot(Tevap, Tcond):
    # Calculate Carnot COP given evaporation and condensation temperatures in Kelvin.
    dT = np.maximum(Tcond - Tevap, 1e-6)
    return Tevap / dT

def evap_effectiveness(UA, m_dot_w):
    # Calculate heat exchanger effectiveness (epsilon) using NTU method for a simple evaporator.
    Cw = max(m_dot_w * CP_WATER, 1e-9)
    NTU = UA / Cw
    return 1.0 - np.exp(-NTU)

def calculate_capacity_constant_dt(Tw, Ta, Pel, eta_carnot, dT_evap, dT_cond, Qev_max=None):
    # Simple model: Constant approach temperatures.
    Tevap = Tw - dT_evap
    Tcond = Ta + dT_cond
    lift  = Tcond - Tevap
    COP = eta_carnot * (Tevap / lift)
    Q = COP * Pel
    return min(Q, Qev_max) if Qev_max is not None else Q
    
def calculate_capacity_ua_implicit(Tw, Ta, Pel, eta_carnot, UA_evap, UA_cond, Qev_max=None):
    # Intermediate model: Linear dT collapse (UA-based) solved implicitly.
    def f(Q):
        Q = max(Q, 0.0)
        Tevap = Tw - (Q / UA_evap)
        Tcond = Ta + ((Q + Pel) / UA_cond)
        lift = max(Tcond - Tevap, 10.0)
        COP = eta_carnot * (Tevap / lift)
        return COP * Pel - Q
    
    Q_lo, Q_hi = 0.0, 20000.0
    f_lo, f_hi = f(Q_lo), f(Q_hi)
    if f_lo * f_hi > 0:
        return np.nan
    Q = float(brentq(f, Q_lo, Q_hi, maxiter=200))
    return min(Q, Qev_max) if Qev_max is not None else Q

def steady_state_detailed(Tw_in, m_dot_w, Ta, Pel, p: ChillerParams, x0=None):
    # Detailed model: Coupled NTU evap + UA cond + COP.
    eps = evap_effectiveness(p.UA_evap, m_dot_w)
    Cw = m_dot_w * CP_WATER
    
    def Qev_hx(Tevap):
        return eps * Cw * max(Tw_in - Tevap, 0.0)
    
    def Qcond_hx(Tcond):
        return p.UA_cond * max(Tcond - Ta, 0.0)
    
    def Qev_cop(Tevap, Tcond):
        lift = max(Tcond - Tevap, p.lift_min)
        COP = p.eta_carnot * (Tevap / lift)
        return COP * Pel
    
    if x0 is None:
        x0 = np.array([max(Tw_in - 5.0, p.Tevap_min), min(Ta + 12.0, p.Tcond_max)], float)
    
    def res(x):
        Tevap, Tcond = x
        Tevap = float(np.clip(Tevap, p.Tevap_min, Tw_in))
        Tcond = float(np.clip(Tcond, Ta, p.Tcond_max))
        
        q_hx = Qev_hx(Tevap)
        q_cop = Qev_cop(Tevap, Tcond)
        q = min(q_hx, q_cop, p.Qev_max)
        
        r1 = q_hx - q_cop
        r2 = Qcond_hx(Tcond) - (q + Pel)
        return np.array([r1, r2], float)
    
    sol = root(res, x0, method="hybr")
    if not sol.success:
        return np.nan, np.nan, np.nan, np.nan, x0
    
    Tevap, Tcond = sol.x
    Tevap = float(np.clip(Tevap, p.Tevap_min, Tw_in))
    Tcond = float(np.clip(Tcond, Ta, p.Tcond_max))
    
    q_hx = eps * Cw * max((Tw_in - Tevap), 0.0)
    lift = max(Tcond - Tevap, p.lift_min)
    COP = p.eta_carnot * (Tevap / lift)
    q_cop = COP * Pel
    Qev = min(q_hx, q_cop, p.Qev_max)
    Tw_out = Tw_in - (Qev / max(Cw, 1e-9))
    return Qev, Tevap, Tcond, Tw_out, np.array([Tevap, Tcond], float)

def main():
    # Main execution block: Calibrates parameters and plots comparison of three chiller models.
    # -------------------
    # Given / assumed (Titan 4000)
    # -------------------
    Pel          = 1.25 * u.kW
    Q_rated      = 3.8 * u.kW
    Tw_rated     = u.C2K(22.0)
    Tamb_rated   = u.C2K(32.0)
    # Assumed water flow (not provided)
    m_dot_w = 0.3  # kg/s 

    Ta   = u.C2K(32.0)  # user assumption for the plot

    

    # Choose rating approaches for calibration (to set UA/eta)
    dT_evap_r = 5.0
    dT_cond_r = 10.0
    Tevap_r = Tw_rated - dT_evap_r
    Tcond_r = Tamb_rated + dT_cond_r

    COP_real_r = Q_rated / Pel
    COP_carnot_r = float(cop_carnot(Tevap_r, Tcond_r))
    eta_cal = COP_real_r / COP_carnot_r

    # Calibrate UA_evap (depends on assumed m_dot_w)
    Cw = m_dot_w * CP_WATER
    eps_needed = Q_rated / (Cw * max(Tw_rated - Tevap_r, 1e-9))
    eps_needed = float(np.clip(eps_needed, 1e-6, 0.999999))
    NTU = -np.log(1.0 - eps_needed)
    UA_ev_cal = NTU * Cw

    # Calibrate UA_cond from rating condenser balance
    Qcond_r = Q_rated + Pel
    UA_co_cal = Qcond_r / max((Tcond_r - Tamb_rated), 1e-9)

    params_detailed = ChillerParams(
        UA_evap=UA_ev_cal,
        UA_cond=UA_co_cal,
        eta_carnot=eta_cal,
        lift_min=10.0,
        Qev_max=4.5 * u.kW,  # soft cap
    )

    # -------------------
    # Generate curves
    # -------------------
    Tw_C_range = np.linspace(5, 35, 301)
    Tw_range = u.C2K(Tw_C_range)
    
    Q_detailed = np.empty_like(Tw_range)
    Q_const_dt = np.empty_like(Tw_range)
    Q_ua_implicit = np.empty_like(Tw_range)

    x0 = None
    for i, t in enumerate(Tw_range):
        # Constant dT model
        Q_const_dt[i] = calculate_capacity_constant_dt(t, Ta, Pel, eta_cal, dT_evap_r, dT_cond_r, params_detailed.Qev_max)
        
        # UA Implicit model
        Q_ua_implicit[i] = calculate_capacity_ua_implicit(t, Ta, Pel, eta_cal, UA_ev_cal, UA_co_cal, params_detailed.Qev_max)
        
        # Detailed model
        Qev, Tev, Tco, Twout, x0 = steady_state_detailed(t, m_dot_w, Ta, Pel, params_detailed, x0=x0)

        Q_detailed[i] = Qev 

    # -------------------
    # Plot
    # -------------------
    plt.figure(figsize=(9, 5.2))
    plt.plot(Tw_C_range, Q_detailed / u.kW, color=mCol(1), linewidth=2, label="Detailed model (NTU evap + UA cond + COP)")
    plt.plot(Tw_C_range, Q_const_dt / u.kW, color=mCol(2), linewidth=2, label="Constant dT model (dT_evap=5K, dT_cond=10K)")
    plt.plot(Tw_C_range, Q_ua_implicit / u.kW, color=mCol(3), linewidth=2, label="UA Implicit model (Linear dT)")

    plt.xlabel("Water temperature (°C)")
    plt.ylabel("Cooling power Q (kW)")
    plt.title("Titan 4000 – Cooling power vs water temperature (operating range 5–35°C)\nAssumed air temperature = 23°C")
    plt.grid(True)
    plt.xlim(5, 35)
    plt.ylim(0, np.nanmax([np.nanmax(Q_detailed / u.kW), np.nanmax(Q_const_dt / u.kW), np.nanmax(Q_ua_implicit / u.kW)]) * 1.05)
    plt.legend()

    plt.text(
        0.02, 0.02,
        f"Calibration: Q=3.8 kW @ Tw=22°C, air=32°C, Pel=1.25 kW\n"
        f"Assumed rating approaches: ΔT_evap=5K, ΔT_cond=10K\n"
        f"Assumed water flow: {m_dot_w*60:.1f} L/min (affects UA_evap)",
        transform=plt.gca().transAxes,
        fontsize=9,
        verticalalignment="bottom"
    )
    plt.show()

    print("Calibration summary (used for all three curves):")
    print(f"  eta_carnot = {eta_cal:.3f}")
    print(f"  UA_evap ≈ {UA_ev_cal:.0f} W/K (depends on assumed water flow)")
    print(f"  UA_cond ≈ {UA_co_cal:.0f} W/K")

if __name__ == "__main__":
    main()

