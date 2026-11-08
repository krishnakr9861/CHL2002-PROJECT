"""
Cubic equations of state for pure real gases.

Generic two-parameter cubic EOS (Smith, Van Ness & Abbott form):

        P = RT/(V - b) - a(T) / [(V + eps*b)(V + sig*b)]

        a(T) = Psi * alpha(Tr, omega) * R^2 Tc^2 / Pc ,   b = Omega * R Tc / Pc

    EOS   sig        eps        Omega     Psi       alpha(Tr)
    vdW   0          0          1/8       27/64     1
    RK    1          0          0.08664   0.42748   Tr^-0.5
    SRK   1          0          0.08664   0.42748   [1 + m(1 - Tr^0.5)]^2, m = 0.480 + 1.574w - 0.176w^2
    PR    1+sqrt2    1-sqrt2    0.07780   0.45724   [1 + k(1 - Tr^0.5)]^2, k = 0.37464 + 1.54226w - 0.26992w^2

Dimensionless groups:   beta = bP/RT ,   q = a/(bRT)

Cubic in Z (obtained by putting V = ZRT/P into the EOS):
    (Z - beta)(Z + eps*beta)(Z + sig*beta) - (Z + eps*beta)(Z + sig*beta) + q*beta*(Z - beta) = 0

Residual (departure) properties at (T, P), with
    I = 1/(sig - eps) * ln[(Z + sig*beta)/(Z + eps*beta)]     (I = beta/Z for vdW, sig = eps)

    H^R / RT = Z - 1 + [dln(alpha)/dln(Tr) - 1] * q * I
    S^R / R  = ln(Z - beta) + [dln(alpha)/dln(Tr)] * q * I
    ln(phi)  = Z - 1 - ln(Z - beta) - q * I

Real-gas properties (reference state: ideal gas at T0 = 298.15 K, P0 = 1 bar, H = S = 0):
    H(T,P) = integral(Cp_ig dT, T0..T) + H^R
    S(T,P) = integral(Cp_ig/T dT, T0..T) - R ln(P/P0) + S^R
"""

import math
import numpy as np

from data import R, get_component

SQRT2 = math.sqrt(2.0)

EOS_PARAMS = {
    #       sig         eps         Omega     Psi
    "vdW": (0.0,        0.0,        1 / 8,    27 / 64),
    "RK":  (1.0,        0.0,        0.08664,  0.42748),
    "SRK": (1.0,        0.0,        0.08664,  0.42748),
    "PR":  (1 + SQRT2,  1 - SQRT2,  0.07780,  0.45724),
}

T0 = 298.15     # K, reference temperature
P0 = 1.0e5      # Pa, reference pressure


# ----------------------------------------------------------------------------
# alpha(Tr) and its logarithmic temperature derivative
# ----------------------------------------------------------------------------
def alpha_function(eos, Tr, omega):
    """Return alpha and dln(alpha)/dln(Tr)."""
    if eos == "vdW":
        return 1.0, 0.0
    if eos == "RK":
        return Tr ** -0.5, -0.5
    if eos == "SRK":
        m = 0.480 + 1.574 * omega - 0.176 * omega ** 2
    elif eos == "PR":
        m = 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2
    else:
        raise ValueError(f"Unknown EOS '{eos}'")
    s = 1.0 + m * (1.0 - math.sqrt(Tr))
    alpha = s * s
    dlna_dlnTr = -m * math.sqrt(Tr) / s
    return alpha, dlna_dlnTr


def eos_constants(comp, T, eos="PR"):
    """a(T) [Pa m6/mol2], b [m3/mol], dln(alpha)/dln(Tr)."""
    sig, eps, Omega, Psi = EOS_PARAMS[eos]
    Tc, Pc, w = comp["Tc"], comp["Pc"], comp["omega"]
    alpha, dlna = alpha_function(eos, T / Tc, w)
    a = Psi * alpha * R ** 2 * Tc ** 2 / Pc
    b = Omega * R * Tc / Pc
    return a, b, dlna


# ----------------------------------------------------------------------------
# Solving the cubic for Z
# ----------------------------------------------------------------------------
def z_roots(beta, q, eos):
    """All physically meaningful real roots (Z > beta) of the cubic, ascending."""
    sig, eps, _, _ = EOS_PARAMS[eos]
    p1 = np.poly1d([1.0, -beta])               # (Z - beta)
    p2 = np.poly1d([1.0, eps * beta])          # (Z + eps*beta)
    p3 = np.poly1d([1.0, sig * beta])          # (Z + sig*beta)
    poly = p1 * p2 * p3 - p2 * p3 + q * beta * p1
    roots = np.roots(poly.coeffs)
    real = sorted(r.real for r in roots if abs(r.imag) < 1e-9 and r.real > beta)
    return real


def solve_Z(T, P, comp, eos="PR", phase="auto"):
    """
    Compressibility factor at (T [K], P [Pa]).

    phase = 'vapor'  -> largest root
            'liquid' -> smallest root
            'auto'   -> if three roots exist, the one with the lowest Gibbs
                        energy (lowest ln phi), i.e. the stable phase.
    Returns (Z, beta, q, dlnalpha, n_roots).
    """
    a, b, dlna = eos_constants(comp, T, eos)
    beta = b * P / (R * T)
    q = a / (b * R * T)
    roots = z_roots(beta, q, eos)
    if not roots:
        raise RuntimeError("No physical root found")
    if phase == "vapor":
        Z = roots[-1]
    elif phase == "liquid":
        Z = roots[0]
    else:
        if len(roots) == 1:
            Z = roots[0]
        else:
            lnphi = [ln_phi_from_Z(z, beta, q, eos) for z in (roots[0], roots[-1])]
            Z = roots[0] if lnphi[0] < lnphi[1] else roots[-1]
    return Z, beta, q, dlna, len(roots)


# ----------------------------------------------------------------------------
# Residual properties and fugacity
# ----------------------------------------------------------------------------
def I_integral(Z, beta, eos):
    sig, eps, _, _ = EOS_PARAMS[eos]
    if sig == eps:                       # van der Waals
        return beta / (Z + eps * beta)
    return math.log((Z + sig * beta) / (Z + eps * beta)) / (sig - eps)


def ln_phi_from_Z(Z, beta, q, eos):
    return Z - 1.0 - math.log(Z - beta) - q * I_integral(Z, beta, eos)


def residual_properties(Z, beta, q, dlna, eos):
    """Return HR/RT, SR/R, ln(phi)."""
    I = I_integral(Z, beta, eos)
    HR_RT = Z - 1.0 + (dlna - 1.0) * q * I
    SR_R = math.log(Z - beta) + dlna * q * I
    lnphi = Z - 1.0 - math.log(Z - beta) - q * I
    return HR_RT, SR_R, lnphi


# ----------------------------------------------------------------------------
# Ideal-gas contributions (SVNA Cp form)
# ----------------------------------------------------------------------------
def ideal_gas_H(comp, T, Tref=T0):
    """H_ig(T) - H_ig(Tref) in J/mol."""
    A, B, C, D = comp["cp"]
    val = (A * (T - Tref) + B / 2 * (T ** 2 - Tref ** 2)
           + C / 3 * (T ** 3 - Tref ** 3) - D * (1 / T - 1 / Tref))
    return R * val


def ideal_gas_S(comp, T, P, Tref=T0, Pref=P0):
    """S_ig(T,P) - S_ig(Tref,Pref) in J/(mol K)."""
    A, B, C, D = comp["cp"]
    val = (A * math.log(T / Tref) + B * (T - Tref) + C / 2 * (T ** 2 - Tref ** 2)
           - D / 2 * (1 / T ** 2 - 1 / Tref ** 2))
    return R * val - R * math.log(P / Pref)


def ideal_gas_Cp(comp, T):
    A, B, C, D = comp["cp"]
    return R * (A + B * T + C * T ** 2 + D / T ** 2)


# ----------------------------------------------------------------------------
# Everything at one state point
# ----------------------------------------------------------------------------
def state(gas, T, P, eos="PR", phase="auto"):
    """
    Full property set for a pure component.

    gas : component name (see data.py) or a component dict
    T   : temperature, K
    P   : pressure, Pa
    """
    comp = get_component(gas) if isinstance(gas, str) else gas
    Z, beta, q, dlna, nroots = solve_Z(T, P, comp, eos, phase)
    HR_RT, SR_R, lnphi = residual_properties(Z, beta, q, dlna, eos)

    HR = HR_RT * R * T
    SR = SR_R * R
    H = ideal_gas_H(comp, T) + HR
    S = ideal_gas_S(comp, T, P) + SR
    phi = math.exp(lnphi)
    a, b, _ = eos_constants(comp, T, eos)

    return {
        "gas": comp["name"], "eos": eos, "T": T, "P": P,
        "Tr": T / comp["Tc"], "Pr": P / comp["Pc"],
        "a": a, "b": b, "A": q * beta, "B": beta,
        "n_roots": nroots, "Z": Z, "V": Z * R * T / P,
        "HR": HR, "SR": SR, "H": H, "S": S,
        "phi": phi, "f": phi * P,
    }


# ----------------------------------------------------------------------------
# Bonus: saturation pressure from equal fugacities (phi_L = phi_V)
# ----------------------------------------------------------------------------
def saturation_pressure(gas, T, eos="PR", tol=1e-10, maxit=500):
    comp = get_component(gas) if isinstance(gas, str) else gas
    Tc, Pc, w = comp["Tc"], comp["Pc"], comp["omega"]
    if T >= Tc:
        raise ValueError("T must be below Tc for a saturation pressure")
    P = Pc * math.exp(5.373 * (1 + w) * (1 - Tc / T))     # Wilson estimate
    for _ in range(maxit):
        a, b, _ = eos_constants(comp, T, eos)
        beta = b * P / (R * T)
        q = a / (b * R * T)
        roots = z_roots(beta, q, eos)
        if len(roots) < 2:
            # outside the 3-root region: nudge P toward it
            P *= 0.9 if roots[0] < 0.3 else 1.1
            continue
        phiL = math.exp(ln_phi_from_Z(roots[0], beta, q, eos))
        phiV = math.exp(ln_phi_from_Z(roots[-1], beta, q, eos))
        Pnew = P * phiL / phiV
        if abs(Pnew / P - 1) < tol:
            return Pnew, roots[0], roots[-1]
        P = Pnew
    raise RuntimeError("Saturation pressure did not converge")
