import math

def cop_carnot(Tevap_K: float, Tcond_K: float) -> float:
    if Tcond_K == Tevap_K:
        return float('inf')
    return Tevap_K / (Tcond_K - Tevap_K)

def calculate_chiller_performance(
    Tamb_C: float, 
    P_kW: float, 
    Q_rated_kW: float, 
    Tamb_rated_C: float, 
    Twater_rated_C: float, 
    dT_evap: float, 
    dT_cond: float
):
    """
    Calculate chiller performance and return eta, plot_data, and reference_points.
    """
    # Calculate eta
    Tevap_rated_K = Twater_rated_C - dT_evap + 273.15
    Tcond_rated_K = Tamb_rated_C + dT_cond + 273.15
    COP_carnot_rated = cop_carnot(Tevap_rated_K, Tcond_rated_K)
    COP_real_rated = Q_rated_kW / P_kW
    eta = COP_real_rated / COP_carnot_rated

    # Data for plot
    plot_data = []
    # Range 8 to 25 step 0.5
    for i in range(0, 35): 
        t = 8.0 + i * 0.5
        if t > 25.0: break
        
        Tevap_K = t - dT_evap + 273.15
        Tcond_K = Tamb_C + dT_cond + 273.15
        cop = cop_carnot(Tevap_K, Tcond_K)
        cop_real = eta * cop
        q = P_kW * cop_real
        
        if not math.isfinite(q) or q <= 0:
            q_val = None
        else:
            q_val = q

        plot_data.append({"Twater_C": round(t, 1), "Q_kW": q_val if q_val is not None else 0.0})

    # Reference points
    reference_points = []
    for t in [8, 11, 15]:
        te = t - dT_evap + 273.15
        tc = Tamb_C + dT_cond + 273.15
        q = P_kW * eta * cop_carnot(te, tc)
        reference_points.append({"Twater_C": t, "Q_kW": q})

    return {
        "eta": eta,
        "plot_data": plot_data,
        "reference_points": reference_points
    }
