"""
Pure-component data used by the program.

Critical constants (Tc, Pc, omega) and ideal-gas heat-capacity coefficients
are taken from Smith, Van Ness & Abbott, "Introduction to Chemical
Engineering Thermodynamics", Appendix B (Table B.1) and Appendix C (Table C.1).

Ideal-gas heat capacity form (SVNA Table C.1):
    Cp_ig / R = A + B*T + C*T**2 + D*T**-2        (T in K)
"""

R = 8.314462618  # J/(mol K)

# name: (Tc [K], Pc [bar], omega, A, B, C, D)
COMPONENTS = {
    "methane":        (190.6,  45.99, 0.012, 1.702, 9.081e-3, -2.164e-6,  0.0),
    "ethane":         (305.3,  48.72, 0.100, 1.131, 19.225e-3, -5.561e-6, 0.0),
    "propane":        (369.8,  42.48, 0.152, 1.213, 28.785e-3, -8.824e-6, 0.0),
    "n-butane":       (425.1,  37.96, 0.200, 1.935, 36.915e-3, -11.402e-6, 0.0),
    "nitrogen":       (126.2,  34.00, 0.038, 3.280, 0.593e-3,  0.0,  0.040e5),
    "oxygen":         (154.6,  50.43, 0.022, 3.639, 0.506e-3,  0.0, -0.227e5),
    "carbon dioxide": (304.2,  73.83, 0.224, 5.457, 1.045e-3,  0.0, -1.157e5),
    "ammonia":        (405.7, 112.80, 0.253, 3.578, 3.020e-3,  0.0, -0.186e5),
    "water":          (647.1, 220.55, 0.345, 3.470, 1.450e-3,  0.0,  0.121e5),
}


def get_component(name):
    key = name.strip().lower()
    if key in ("co2",):
        key = "carbon dioxide"
    if key in ("ch4",):
        key = "methane"
    if key in ("n2",):
        key = "nitrogen"
    if key in ("h2o", "steam"):
        key = "water"
    if key not in COMPONENTS:
        raise KeyError(f"Unknown component '{name}'. Available: {', '.join(COMPONENTS)}")
    Tc, Pc_bar, w, A, B, C, D = COMPONENTS[key]
    return {
        "name": key,
        "Tc": Tc,
        "Pc": Pc_bar * 1e5,  # Pa
        "omega": w,
        "cp": (A, B, C, D),
    }
