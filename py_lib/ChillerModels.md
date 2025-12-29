# Chiller Modeling Methods

This document describes the three chiller capacity prediction methods implemented in `Chiller.py`.

## Overview

All three models predict the cooling capacity of a vapor-compression chiller given operating conditions. They differ in complexity and accuracy, representing trade-offs between computational simplicity and physical realism.

### Common Inputs
- **Tw**: Water temperature entering the evaporator (K)
- **Ta**: Ambient air temperature at the condenser (K)
- **Pel**: Electrical power input to the compressor (W)
- **eta_carnot**: Efficiency factor relative to ideal Carnot COP (0-1)
- **Qev_max**: Maximum equipment cooling capacity (W)

### Common Physics
All models use the modified Carnot relationship:

```
COP_real = eta_carnot × COP_carnot
COP_carnot = T_evap / (T_cond - T_evap)
```

Where:
- `T_evap`: Refrigerant evaporation temperature (K)
- `T_cond`: Refrigerant condensation temperature (K)

---

## Method 1: Constant ΔT Model

### Function
`calculate_capacity_constant_dt(Tw, Ta, Pel, eta_carnot, dT_evap, dT_cond, Qev_max=None)`

### Description
The simplest model assumes fixed temperature approaches at both heat exchangers.

### Approach Assumptions
- Evaporator approach: `T_evap = Tw - dT_evap` (constant)
- Condenser approach: `T_cond = Ta + dT_cond` (constant)

### Algorithm
1. Calculate refrigerant temperatures using constant approaches
2. Calculate temperature lift: `lift = T_cond - T_evap`
3. Calculate COP: `COP = eta_carnot × (T_evap / lift)`
4. Calculate capacity: `Q = COP × Pel`
5. Apply equipment limit: `Q = min(Q, Qev_max)`

### Advantages
- Extremely fast (direct calculation)
- No iteration required
- Useful for quick estimates

### Limitations
- Ignores heat exchanger physics
- Approach temperatures don't vary with load
- Less accurate at off-design conditions
- Requires empirical selection of dT_evap and dT_cond

### Typical Parameters
- dT_evap ≈ 5 K
- dT_cond ≈ 10 K

---

## Method 2: UA Implicit Model

### Function
`calculate_capacity_ua_implicit(Tw, Ta, Pel, eta_carnot, UA_evap, UA_cond, Qev_max=None)`

### Description
An intermediate model where approach temperatures vary linearly with heat transfer rate based on UA values.

### Heat Transfer Relations
- Evaporator: `Q = UA_evap × (Tw - T_evap)`
- Condenser: `Q + Pel = UA_cond × (T_cond - Ta)`

### Algorithm
Solves the coupled system implicitly using root finding:

1. Define residual function f(Q):
   - From evaporator: `T_evap = Tw - Q/UA_evap`
   - From condenser: `T_cond = Ta + (Q + Pel)/UA_cond`
   - Calculate lift (with minimum enforced)
   - Calculate COP
   - Residual: `f(Q) = COP × Pel - Q`

2. Find root using Brent's method in range [0, 20000] W
3. Apply equipment limit: `Q = min(Q, Qev_max)`

### Advantages
- Accounts for heat exchanger performance
- Approach temperatures vary with capacity
- More realistic than constant ΔT
- Still relatively fast (1D root finding)

### Limitations
- Assumes linear ΔT relationships (not exact for real heat exchangers)
- Doesn't account for water flow rate or effectiveness
- No energy balance on water side

### Typical Parameters
- UA_evap ≈ 1000-5000 W/K
- UA_cond ≈ 300-1000 W/K

---

## Method 3: Detailed Model (NTU Method)

### Function
`steady_state_detailed(Tw_in, m_dot_w, Ta, Pel, p: ChillerParams, x0=None)`

### Description
The most sophisticated model using NTU-effectiveness method for the evaporator and full energy balances.

### Physical Sub-Models

#### Evaporator (NTU Method)
- Heat capacity rate: `Cw = m_dot_w × Cp_water`
- NTU: `NTU = UA_evap / Cw`
- Effectiveness: `ε = 1 - exp(-NTU)`
- Capacity: `Q_hx = ε × Cw × (Tw_in - T_evap)`

#### Condenser (UA Method)
- Heat rejection: `Q_cond = UA_cond × (T_cond - Ta)`
- Energy balance: `Q_cond = Q_evap + Pel`

#### Compressor
- COP: `COP = eta_carnot × T_evap / (T_cond - T_evap)`
- Capacity: `Q_cop = COP × Pel`

### Algorithm
Solves for [T_evap, T_cond] using 2D root finding:

1. **Residual 1**: Match evaporator HX and COP capacities
   - `r1 = Q_hx - Q_cop`
   
2. **Residual 2**: Condenser energy balance
   - `r2 = Q_cond_hx - (Q + Pel)`

3. **Limiting**: Final capacity is minimum of:
   - Heat exchanger capacity (Q_hx)
   - Thermodynamic capacity (Q_cop)
   - Equipment maximum (Qev_max)

4. **Output water temperature**:
   - `Tw_out = Tw_in - Q / (m_dot_w × Cp_water)`

### Safety Constraints
- `T_evap_min`: Freeze protection (typically 0.5°C)
- `T_cond_max`: Safety limit (typically 70°C)
- `lift_min`: Minimum temperature lift for operation (typically 10 K)

### Advantages
- Most physically accurate
- Accounts for water flow rate effects
- Proper NTU-effectiveness for evaporator
- Full energy balance on both sides
- Returns water outlet temperature
- Respects equipment operating limits

### Limitations
- Requires iteration (2D root finding)
- Needs more parameters (UA_evap, UA_cond, m_dot_w)
- Slightly slower computation
- May fail to converge in extreme conditions

### Typical Parameters
- m_dot_w ≈ 0.5-2.0 kg/s (typical water flow)
- UA_evap ≈ 1000-5000 W/K
- UA_cond ≈ 300-1000 W/K
- eta_carnot ≈ 0.3-0.5
- lift_min = 10 K

---

## Important Note on Flow Rate Sensitivity

**Method 2 (UA Implicit) does not include water flow rate in its physics—only in the calibration of UA_evap.**

- When calibrating at a low water flow rate, UA_evap becomes small (since it is based on NTU-effectiveness).
- However, Method 2 assumes a simple linear relationship between Q and (Tw - Tevap), and does not account for the diminishing effectiveness at low flow.
- As a result, Method 2 will not match the calibration point at low flow rates, while Method 3 (NTU) will.
- This is a limitation of the linear UA model: it cannot capture the true drop in performance at low flow, because it lacks the NTU-effectiveness physics.

**Recommendation:** For accurate results at low flow rates, use Method 3 (NTU-effectiveness). Method 2 is best suited for moderate-to-high flow rates where effectiveness is close to 1.

---

## Parameter Calibration

All models can be calibrated from manufacturer's rating data:

### Given Rating Point Data
- Q_rated: Rated cooling capacity
- Tw_rated: Rated water inlet temperature
- Ta_rated: Rated ambient temperature
- Pel: Electrical power consumption

### Calibration Steps

1. **Choose rating approaches** (for UA models):
   ```python
   dT_evap_r = 5.0  # K
   dT_cond_r = 10.0  # K
   ```

2. **Calculate rating point temperatures**:
   ```python
   T_evap_r = Tw_rated - dT_evap_r
   T_cond_r = Ta_rated + dT_cond_r
   ```

3. **Calculate eta_carnot** (same for all models):
   ```python
   COP_real = Q_rated / Pel
   COP_carnot = T_evap_r / (T_cond_r - T_evap_r)
   eta_carnot = COP_real / COP_carnot
   ```

4. **For UA models, calibrate heat exchangers**:

   #### Finding UA_evap (Evaporator)
   
   ⚠️ **Two approaches exist, but current implementation uses NTU method for both:**
   
   **Approach A: Simple Linear (could be used for Method 2):**
   
   ```python
   # From: Q = UA_evap × (Tw - T_evap)
   # Solving: UA_evap = Q / ΔT
   UA_evap = Q_rated / (Tw_rated - T_evap_r)
   ```
   - Independent of water flow rate
   - Simpler but less physically accurate
   
   **Approach B: NTU Method (current implementation for both methods):**
   
   ```python
   # Water-side heat capacity rate
   Cw = m_dot_w × Cp_water  # Requires assumed water flow rate
   
   # Required effectiveness at rating point
   # From: Q = ε × Cw × (Tw_in - T_evap)
   ε_needed = Q_rated / (Cw × (Tw_rated - T_evap_r))
   
   # Convert effectiveness to NTU
   # For counter-flow HX: ε = 1 - exp(-NTU)
   # Solving: NTU = -ln(1 - ε)
   NTU = -log(1 - ε_needed)
   
   # Calculate UA from definition: NTU = UA / Cw
   UA_evap = NTU × Cw
   ```
   - More physically accurate
   - **Depends on assumed water flow rate**
   
   **Current Implementation**: Both Method 2 and Method 3 use Approach B (NTU calibration). This means:
   - Changing the assumed flow rate during calibration changes UA_evap
   - This affects predictions from **both** Method 2 and Method 3
   - Method 2 doesn't explicitly use flow rate in its physics, but inherits the flow-dependent UA_evap from calibration
   
   #### Finding UA_cond (Condenser)
   
   **Same for both Method 2 and Method 3:**
   
   Both use simple UA model with energy balance:
   
   ```python
   # Total heat rejection at rating point
   Q_cond_rated = Q_rated + Pel
   
   # From: Q_cond = UA_cond × (T_cond - Ta)
   # Solving: UA_cond = Q_cond / ΔT
   UA_cond = Q_cond_rated / (T_cond_r - Ta_rated)
   ```
   
   **Physical meaning**: 
   - UA_evap represents the overall heat transfer coefficient × area for the water-refrigerant evaporator
   - UA_cond represents the overall heat transfer coefficient × area for the air-refrigerant condenser
   - Both are in units of W/K (heat transfer rate per degree temperature difference)

---

## Model Selection Guide

| Use Case | Recommended Model |
|----------|-------------------|
| Quick screening calculations | Constant ΔT |
| Preliminary design | UA Implicit |
| Detailed system simulation | Detailed (NTU) |
| Real-time control | Constant ΔT or UA Implicit |
| High accuracy required | Detailed (NTU) |
| Unknown water flow rate | Constant ΔT or UA Implicit |
| Part-load performance | Detailed (NTU) |

---

## Example Usage

```python
from Chiller import ChillerParams, calculate_capacity_constant_dt, \
                    calculate_capacity_ua_implicit, steady_state_detailed
from matLib import u

# Define operating conditions
Tw = u.C2K(20.0)  # 20°C water inlet
Ta = u.C2K(30.0)  # 30°C ambient
Pel = 1.25e3      # 1.25 kW electrical input
m_dot_w = 1.0     # 1 kg/s water flow

# Calibrated parameters
eta_carnot = 0.45
UA_evap = 3000.0  # W/K
UA_cond = 600.0   # W/K
Qev_max = 4500.0  # W

# Method 1: Constant ΔT
Q1 = calculate_capacity_constant_dt(Tw, Ta, Pel, eta_carnot, 
                                     dT_evap=5.0, dT_cond=10.0, 
                                     Qev_max=Qev_max)

# Method 2: UA Implicit
Q2 = calculate_capacity_ua_implicit(Tw, Ta, Pel, eta_carnot, 
                                     UA_evap, UA_cond, 
                                     Qev_max=Qev_max)

# Method 3: Detailed
params = ChillerParams(
    UA_evap=UA_evap,
    UA_cond=UA_cond,
    eta_carnot=eta_carnot,
    Qev_max=Qev_max
)
Q3, Tevap, Tcond, Tw_out, _ = steady_state_detailed(
    Tw, m_dot_w, Ta, Pel, params
)

print(f"Constant ΔT: {Q1/1e3:.2f} kW")
print(f"UA Implicit: {Q2/1e3:.2f} kW")
print(f"Detailed:    {Q3/1e3:.2f} kW (Tw_out = {u.K2C(Tw_out):.1f}°C)")
```

---

## References

- NTU-effectiveness method: Incropera & DeWitt, "Fundamentals of Heat and Mass Transfer"
- Vapor compression cycles: ASHRAE Handbook - Refrigeration
- Heat exchanger design: Kays & London, "Compact Heat Exchangers"
