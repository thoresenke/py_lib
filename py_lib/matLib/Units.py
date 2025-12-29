import math


# Placeholder for uBase if you want inheritance, 
# otherwise we define base units directly in u.
class uBase:
    pass

# Placeholder for UnitType since the definition wasn't provided.
# This allows the code to compile and run.
class UnitType:
    def __init__(self, value, string_rep):
        self.value = value
        self.Str = string_rep

    # Allow division (val / UnitType) to work like (val / UnitType.value)
    def __rtruediv__(self, other):
        return other / self.value

    # Mock enums for the class definition
    Pressure = None # Defined later after u is built
    ECD = None
    Flowrate = None
    Force = None
    Torque = None
    Length = None
    Diameter = None
    Ratio = None
    Volume = None
    Velocity = None
    Angle = None
    Density = None
    NoneType = None
    Rotation = None
    WOB = None
    ROP = None
    Area = None
    Nozzle = None
    Temperature = None
    TemperatureGradient = None
    Time = None
    Stress = None
    ThermalConductivity = None
    SpecHeatCapacity = None
    Fann = None
    Viscosity = None
    DeltaPressure = None
    MassFlowRateTon = None
    PowerMega = None
    Hookload = None
    Mass = None

class Units(uBase):
    """
    Description: unit handling class
    python port of KalTek MATLAB class
    """

    # ============ SI BASE UNITS ============
    # We define these to 1.0 so the logic below works
    m   = 1.0
    kg  = 1.0
    s   = 1.0
    A   = 1.0
    K   = 1.0
    mol = 1.0
    mmol = 1e-3 * mol
    cd  = 1.0

    # ============ DERIVED UNITS ============
    
    # Expressed in terms of other SI units
    # Note: In Python class body, we refer to 'm', not 'u.m'
    rad  = m / m           # Angle Radian
    sr   = m**2 / m**2     # Solid angle Steradian
    Hz   = 1.0 / s         # Frequency hertz
    N    = kg * m / (s**2) # Force newton
    Pa   = N / m**2        # Pressure pascal
    J    = N * m           # Energy joule
    W    = J / s           # Power watt
    C    = s * A           # Electric charge coulomb
    V    = W / A           # Voltage volt
    F    = C / V           # Electric capacitance farad
    ohm  = V / A           # Electric resistance ohm
    S    = A / V           # Electrical conductance siemens
    Wb   = V * s           # Magnetic flux weber
    T    = Wb / m**2       # Magnetic field strength tesla
    H    = Wb / A          # Inductance henry
    lm   = cd * sr         # Luminous flux lumen
    lx   = lm / m**2       # Illuminance lux
    Bq   = 1.0 / s         # Radioactivity becquerel
    Gy   = J / kg          # Absorbed dose gray
    Sv   = J / kg          # Equivalent dose sievert
    kat  = mol / s         # Catalytic activity katal
    
    # ------- length ----
    km  = 1e3 * m
    cm  = 1e-2 * m
    dm  = 1e-1 * m
    mm  = 1e-3 * m
    um  = 1e-6 * m
    nm  = 1e-9 * m
    ang = 1e-10 * m
    inch = 2.54 * cm       # 'in' is a reserved keyword in Python, renamed to 'inch'
    mil = 1e-3 * inch
    ft  = 12 * inch
    yd  = 3 * ft
    mi  = 5280 * ft
    a0  = 0.529e-10 * m
    
    # ------- Volume -------
    cc    = (cm)**3
    L     = 1000 * cc
    mL    = cc
    gal   = 231 * inch**3  # US gallon
    quart = gal / 4
    floz  = gal / 128
    pint  = 16 * floz      # pint (USA, wet)
    bbl   = 42 * gal       # US oil barrel
    
    # ----- mass ---------
    gm    = 1e-3 * kg      # gram
    mg    = 1e-3 * gm      # milligram
    lb    = 0.45359237 * kg # pound mass
    lbm   = 0.45359237 * kg # pound mass
    oz    = (1.0/16.0) * lb
    amu   = 1.66e-27 * kg
    tonne = 1000 * kg      # metric ton
    tn    = 2000 * lbm     # Short ton
    
    # densities
    ppg = lb / gal         # pound per gallon
    sg  = 1000 * kg / m**3 # specific gravity       

    # ---- time -------
    ms  = 1e-3 * s
    us  = 1e-6 * s
    ns  = 1e-9 * s
    ps  = 1e-12 * s
    min = 60 * s
    hr  = 60 * min
    day = 24 * hr
    yr  = 365.242199 * day
    
    # ---- angular units --------------
    deg = 2 * math.pi / 360 # deg -> rad
    rev = 2 * math.pi
    
    # ---- frequency ----
    kHz = 1e3 * Hz
    MHz = 1e6 * Hz
    GHz = 1e9 * Hz
    rpm = rev / min         # rpm -> rad/s
    spm = rpm               # strokes per minute
    
    # ---- force -------
    dyne = 1e-5 * N
    lbf  = 4.44822 * N
    klbf = 1000 * lbf
    tf   = 1000 * kg * 9.81 # Metric ton - force (approx gravity)
    tnf  = 2000 * lbf       # Short ton - force
    kN   = 1000 * N 
    
    # ---- Speed -------
    kmh = km / hr
    
    # ---- torque -------
    Nm  = N * m
    kNm = 1000 * Nm
    
    # ---- flow ---------
    lpm = L / min           # lpm -> m3/s
    gpm = gal / min         # gpm -> m3/s
    m3phr = m**3 / hr
    # ---- fluid viscosity data ---------
    fann = 1.065 * lbf / (100 * ft**2) # from fann manual
    g300 = 511              # shear rate (1/s) at fann speed = 300 rpm
    
    kGam = g300 / 300.0     # shear rate coefficient (1/s)/rpm
    cP   = 1e-3 * Pa * s    # centiPose
    cSt  = 1 * mm**2 / s    # centiStokes
    
    # ----- energy -----
    MJ   = 1e6 * J
    kJ   = 1e3 * J
    mJ   = 1e-3 * J
    uJ   = 1e-6 * J
    nJ   = 1e-9 * J
    eV   = 1.6022e-19 * J
    BTU  = 1.0550559e3 * J
    kWh  = 3.6e6 * J
    cal  = 4.1868 * J
    kCal = 1e3 * cal
    
    # ---- temperature ---
    mK = 1e-3 * K
    uK = 1e-6 * K
    nK = 1e-9 * K
    
    # %%-----thermal
    k_cond = W / (m * K)

    # Note: Relative offset temperatures (like 0 degC) are tricky 
    # in multiplication classes. Usually represented as 1 unit interval.
    # We set them to 1.0 (Interval) here.
    degF = 1.0 
    degC = 1.0
    
    # ---- pressure -----
    torr  = 133.322 * Pa
    mtorr = 1e-3 * torr
    bar   = 1e5 * Pa
    mbar  = 1e-3 * bar
    atm   = 1.013e5 * Pa
    psi   = 6.8947548e3 * Pa
    
    # ----- power --- ---
    MW = 1e6 * W
    kW = 1e3 * W
    mW = 1e-3 * W
    uW = 1e-6 * W
    nW = 1e-9 * W
    pW = 1e-12 * W
    hp = 745.69987 * W
    
    # ------ Voltage -----
    kV = 1e3 * V
    mV = 1e-3 * V
    uV = 1e-6 * V
    
    # ----- Current ------
    mA = 1e-3 * A
    uA = 1e-6 * A
    nA = 1e-9 * A
    
    # ----magnetic field -----
    gauss = 1e-4 * T
    
    # %% Choke data
    rho_ref = 998.203 * kg / m**3  # density of water @20 degC
    # CvUS calculation (Using math.sqrt)
    CvUS = math.sqrt(rho_ref * gal / lb) * gal / min * math.sqrt(lb / gal / psi)
    
    # %% unit handle 
    # (Using integers/floats here to act as keys if UnitType wasn't defined, 
    # but strictly speaking these should match your UnitType enum)
    # To make the code below valid, we assign these *after* defining the values
    # In a real python app, you might define these in the __init__ or separate file.
    pass

    @staticmethod
    def lpmd():
        return 51

    @staticmethod
    def C2K(C):
        """ Kelvin=C2K(C) """
        try:
            import numpy as np
            return np.asarray(C) + 273.15
        except (ImportError, TypeError):
            try:
                return [c + 273.15 for c in C]
            except TypeError:
                return C + 273.15

    @staticmethod
    def K2C(K):
        """ Celsius=K2C(K) """
        try:
            import numpy as np
            return np.asarray(K) - 273.15
        except (ImportError, TypeError):
            try:
                return [k - 273.15 for k in K]
            except TypeError:
                return K - 273.15

    @staticmethod
    def F2C(F):
        """ Celsius=F2C(F) """
        try:
            import numpy as np
            return (np.asarray(F) - 32.0) / 1.8
        except (ImportError, TypeError):
            try:
                return [(f - 32.0) / 1.8 for f in F]
            except TypeError:
                return (F - 32.0) / 1.8

    @staticmethod
    def C2F(C):
        """ Fahrenheit=C2F(C) """
        try:
            import numpy as np
            return np.asarray(C) * 1.8 + 32.0
        except (ImportError, TypeError):
            try:
                return [c * 1.8 + 32.0 for c in C]
            except TypeError:
                return C * 1.8 + 32.0

    @staticmethod
    def str(val, unit_type, inspec=None):
        """
        str=prt(val,UnitType,spec)
        automatic create string of number with units
        """
        spec = inspec
        
        if spec is None:
            # Mimic switch/case with if/elif
            # Note: We are comparing against the objects defined below
            if unit_type == UnitType.Length:
                spec = '{:.0f}'
            elif unit_type == UnitType.Fann:
                spec = '{:.1f}'
            elif unit_type == UnitType.Velocity:
                spec = '{:.2f}'
            else:
                spec = ' {:.1f}'
        
        # Python string formatting
        # Assuming val is a list/array or single number. 
        # For simplicity, handling single number here.
        try:
            # If iterable (list/tuple)
            formatted_vals = [spec.format(v / unit_type.value) for v in val]
            return ",".join(formatted_vals) + " " + unit_type.Str
        except TypeError:
            # Single number
            return (spec.format(val / unit_type.value)) + " " + unit_type.Str

    @staticmethod
    def st(val, unit_type, inspec=None):
        """
        In Python, getting the 'inputname' (variable name) of the caller 
        is not standard/reliable. Returning just the value string.
        """
        if inspec:
            s = Units.str(val, unit_type, inspec)
        else:
            s = Units.str(val, unit_type)
        
        return s

    @staticmethod
    def setSymBase():
        # Placeholder for file operations
        print("Setting SymBase - Python logic would go here (shutil.copy)")

    @staticmethod
    def setNumBase():
        # Placeholder for file operations
        print("Setting NumBase - Python logic would go here")

# ==========================================
# Initialize UnitType Helpers (Post-Class)
# ==========================================
# We define these here because they rely on 'u' values
UnitType.Pressure = UnitType(Units.Pa, "Pa")
UnitType.ECD = UnitType(Units.sg, "sg") # Example
UnitType.Flowrate = UnitType(Units.lpm, "lpm")
UnitType.Force = UnitType(Units.kN, "kN")
UnitType.Torque = UnitType(Units.kNm, "kNm")
UnitType.Length = UnitType(Units.m, "m")
UnitType.Diameter = UnitType(Units.inch, "in")
UnitType.Ratio = UnitType(1, "")
UnitType.Volume = UnitType(Units.L, "L")
UnitType.Velocity = UnitType(Units.kmh, "km/h")
UnitType.Angle = UnitType(Units.deg, "deg")
UnitType.Density = UnitType(Units.sg, "sg")
UnitType.Rotation = UnitType(Units.rpm, "rpm")
UnitType.Fann = UnitType(Units.fann, "fann")
# Add others as needed...

def test_units():
    print("--- Testing KalTek Unit Class ---")
    
    # Test 1: Basic Length Calculation
    # MATLAB Equivalent: A = 4 * u.km
    A = 4 * Units.km
    print(f"4 km in base units (m): {A}")
    assert A == 4000.0, "Basic km to m conversion failed"
    
    # Test 2: Inverse Calculation
    # MATLAB Equivalent: B = A / u.km
    B = A / Units.km
    print(f"4000m converted back to km: {B}")
    assert B == 4.0, "Inverse conversion failed"

    # Test 3: String Formatting using u.str
    # This tests if UnitType integration works with the static string method
    formatted_len = Units.str(A, UnitType.Length)
    print(f"Formatted String (Length): {formatted_len}")
    
    # Test 4: Viscosity/Fann reading formatting
    # Fann unit is 1.065*lbf/(100ft^2)
    fann_val = 10 * Units.fann
    formatted_fann = Units.str(fann_val, UnitType.Fann)
    print(f"Formatted Fann (10 units): {formatted_fann}")
    
    # Test 5: Temperature
    # Boiling point of water check
    f_boil = 212
    c_boil = Units.F2C(f_boil)
    print(f"212 F in Celsius: {c_boil}")
    assert abs(c_boil - 100.0) < 0.001, "Temp conversion inaccurate"

    # Test 6: Flowrate
    # Boiling point of water check
    Q =  500*Units.lpm

    print(f"500 lpm in m^3/s: {Q}")
    assert abs(Q - 0.00833333333333333) < 0.001, "Flowrate  inaccurate"



    print("--- All Tests Passed ---")

   # p=500*uu.Pressure
   # print(p)


if __name__ == "__main__":
    test_units()
