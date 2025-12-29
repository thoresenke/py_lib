# CO₂ Degassing Theory

This document outlines the theoretical models and equations used in the `degassing_si.py` module for simulating CO₂ removal from water.

## 1. Equilibrium CO₂ Concentration ($C^*$)

```python
def calculate_equilibrium_co2(temp_k: float, salinity_ppt: float, pressure_pa: float = 101325.0) -> float:
```


The equilibrium concentration of CO₂ in water, denoted as $C^*$, is determined by Henry's Law, which states that the solubility of a gas in a liquid is directly proportional to the partial pressure of that gas above the liquid.

$$ C^* = K_H \cdot pCO_2 $$

Where:
*   $C^*$ is the equilibrium concentration (e.g., mg/L or mol/kg).
*   $K_H$ is Henry's law constant, which depends on temperature and salinity.
*   $pCO_2$ is the partial pressure of CO₂ in the atmosphere (approx. 420 ppm).

### Simplified Engineering Approximation
In the code, we use a simplified empirical fit based on the Weiss (1974) equation to estimate $C^*$ as a function of Temperature ($T$) and Salinity ($S$).

$$ C^* \approx C_{ref} \cdot e^{-0.03(T - 15)} \cdot (1 - 0.005(S - 34)) $$

Where:
*   $C_{ref} \approx 0.6$ mg/L (Base solubility at 15°C, 34 ppt, and atmospheric $pCO_2$).
*   $T$ is temperature in °C.
*   $S$ is salinity in ppt.

> **Note:** This is a simplified engineering approximation. For high-precision scientific applications, the full Weiss (1974) equation should be implemented.

## 2. Carbonate System and pH

The pH of the water is estimated based on the concentration of dissolved CO₂ and the Total Alkalinity ($TA$). The carbonate system in seawater is governed by the following equilibrium:

$$ CO_2 + H_2O \leftrightarrow H_2CO_3 \leftrightarrow H^+ + HCO_3^- \leftrightarrow 2H^+ + CO_3^{2-} $$

Assuming that at typical aquaculture pH levels (pH < 8.5), the majority of alkalinity is contributed by bicarbonate ions ($HCO_3^-$), we can approximate the relationship using the first dissociation constant ($K_1$):

$$ K_1 = \frac{[H^+][HCO_3^-]}{[CO_2]} $$

Rearranging for $[H^+]$:

$$ [H^+] = K_1 \cdot \frac{[CO_2]}{[HCO_3^-]} $$

Taking the negative logarithm gives the Henderson-Hasselbalch equation form:

$$ pH = pK_1 + \log_{10}\left(\frac{[HCO_3^-]}{[CO_2]}\right) $$

### Implementation Details
*   **$pK_1$ Approximation:** We use a linear fit for $pK_1$ based on Millero et al. (seawater scale):
    $$ pK_1 \approx 6.0 - 0.01(T - 15) $$
*   **Alkalinity:** We assume $[HCO_3^-] \approx TA$ (Total Alkalinity).

## 3. Mass Transfer Model (Degassing)

The removal of CO₂ is modeled using the standard aeration equation for a Continuous Stirred-Tank Reactor (CSTR) or a Plug Flow Reactor (PFR) approximation, often referred to as the "efficiency" or "G/L ratio" model.

The fundamental mass transfer equation is:

$$ \frac{dC}{dt} = K_L a \cdot (C^* - C) $$

For a flow-through system, the outlet concentration $C_{out}$ is calculated as:

$$ C_{out} = C^* + (C_{in} - C^*) \cdot e^{-K_L a \cdot \tau} $$

Where:
*   $C_{in}$ is the inlet CO₂ concentration.
*   $\tau$ is the retention time ($V/Q_{water}$).
*   $K_L a$ is the volumetric mass transfer coefficient.

### Efficiency Factor ($k$)
In many engineering applications, $K_L a$ is assumed to be proportional to the specific airflow rate ($Q_{air}/V$). This leads to the Gas-to-Liquid ($G/L$) ratio model:

$$ K_L a \cdot \tau \approx k \cdot \frac{Q_{air}}{V} \cdot \frac{V}{Q_{water}} = k \cdot \frac{Q_{air}}{Q_{water}} $$

Thus, the final equation used in the simulation is:

$$ C_{out} = C^* + (C_{in} - C^*) \cdot e^{-k \cdot (G/L)} $$

Where:
*   $k$ is an empirical efficiency factor (dimensionless in this form) representing the performance of the specific aerator/degasser unit.
*   $G/L = Q_{air} / Q_{water}$ is the gas-to-liquid ratio.

## References
1.  **Weiss, R. F. (1974).** Carbon dioxide in water and seawater: the solubility of a non-ideal gas. *Marine Chemistry*, 2(3), 203-215.
2.  **Millero, F. J., et al.** (Various publications on the thermodynamics of the carbonate system in seawater).
3.  **Colt, J. (2012).** *Computation of Dissolved Gas Concentration in Water as Functions of Temperature, Salinity, and Pressure*. Elsevier.
