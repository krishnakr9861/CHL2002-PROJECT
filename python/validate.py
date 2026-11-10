"""
Validation of the cubic-EOS program against reference data.

Reference values come from CoolProp (pip install CoolProp), which implements the
same high-accuracy multiparameter reference equations used by the NIST Chemistry
WebBook / REFPROP. Departure functions are compared, so the different reference
states of the two programs do not matter:

    H^R = H(T,P) - H(T, P->0)
    S^R = S(T,P) - S(T, P->0) + R ln(P/P->0)

Writes results/validation.csv and results/psat.csv and prints summary tables.
"""

import csv
import math
import os

import CoolProp.CoolProp as CP
from CoolProp import AbstractState

from data import R
from eos import EOS_PARAMS, state, saturation_pressure

COOLPROP_NAME = {
    "methane": "Methane", "carbon dioxide": "CarbonDioxide", "nitrogen": "Nitrogen",
    "propane": "Propane", "water": "Water", "ethane": "Ethane", "ammonia": "Ammonia",
}

# (gas, T [K], P [bar]) - single-phase gas / supercritical points
CASES = [
    ("methane", 250.0, 50.0),
    ("methane", 300.0, 100.0),
    ("nitrogen", 200.0, 100.0),
    ("nitrogen", 300.0, 200.0),
    ("carbon dioxide", 320.0, 50.0),
    ("carbon dioxide", 350.0, 100.0),
    ("propane", 400.0, 30.0),
    ("water", 573.15, 50.0),
]

PSAT_CASES = [("propane", 300.0), ("carbon dioxide", 280.0), ("water", 373.15),
              ("water", 473.15), ("methane", 150.0), ("ammonia", 300.0)]

P_LOW = 10.0  # Pa, effectively ideal gas


def reference(gas, T, P):
    st = AbstractState("HEOS", COOLPROP_NAME[gas])
    st.update(CP.PT_INPUTS, P, T)
    Z, h, s, phi = st.compressibility_factor(), st.hmolar(), st.smolar(), st.fugacity_coefficient(0)
    st.update(CP.PT_INPUTS, P_LOW, T)
    HR = h - st.hmolar()
    SR = s - st.smolar() + R * math.log(P / P_LOW)
    return {"Z": Z, "HR": HR, "SR": SR, "phi": phi}


def pct(a, b):
    return 100.0 * (a - b) / b


def run(verbose=True):
    os.makedirs("results", exist_ok=True)
    rows = []
    for gas, T, Pbar in CASES:
        ref = reference(gas, T, Pbar * 1e5)
        for eos in EOS_PARAMS:
            s = state(gas, T, Pbar * 1e5, eos)
            rows.append({
                "gas": gas, "T_K": T, "P_bar": Pbar, "eos": eos,
                "Z": s["Z"], "Z_ref": ref["Z"],
                "HR": s["HR"], "HR_ref": ref["HR"],
                "SR": s["SR"], "SR_ref": ref["SR"],
                "phi": s["phi"], "phi_ref": ref["phi"],
                "errZ_pct": pct(s["Z"], ref["Z"]),
                "errHR_pct": pct(s["HR"], ref["HR"]),
                "errSR_pct": pct(s["SR"], ref["SR"]),
                "errphi_pct": pct(s["phi"], ref["phi"]),
            })
    with open("results/validation.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    psat_rows = []
    for gas, T in PSAT_CASES:
        Pref = CP.PropsSI("P", "T", T, "Q", 0, COOLPROP_NAME[gas])
        row = {"gas": gas, "T_K": T, "Psat_ref_bar": Pref / 1e5}
        for eos in EOS_PARAMS:
            row[f"Psat_{eos}_bar"] = saturation_pressure(gas, T, eos)[0] / 1e5
        psat_rows.append(row)
    with open("results/psat.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(psat_rows[0]))
        w.writeheader()
        w.writerows(psat_rows)

    if verbose:
        print("\nPeng-Robinson vs reference (CoolProp / NIST-quality EOS)")
        print(f"{'gas':15s}{'T':>7s}{'P':>6s} | {'Z':>7s}{'Zref':>7s} | {'HR':>8s}{'HRref':>8s}"
              f" | {'SR':>7s}{'SRref':>7s} | {'phi':>6s}{'phiref':>7s}")
        for r in rows:
            if r["eos"] != "PR":
                continue
            print(f"{r['gas']:15s}{r['T_K']:7.1f}{r['P_bar']:6.0f} | {r['Z']:7.4f}{r['Z_ref']:7.4f} | "
                  f"{r['HR']:8.0f}{r['HR_ref']:8.0f} | {r['SR']:7.2f}{r['SR_ref']:7.2f} | "
                  f"{r['phi']:6.4f}{r['phi_ref']:7.4f}")

        print("\nMean absolute % error over all validation points")
        print(f"{'EOS':5s}{'Z':>8s}{'H_R':>8s}{'S_R':>8s}{'phi':>8s}")
        for eos in EOS_PARAMS:
            sub = [r for r in rows if r["eos"] == eos]
            m = lambda k: sum(abs(r[k]) for r in sub) / len(sub)
            print(f"{eos:5s}{m('errZ_pct'):8.2f}{m('errHR_pct'):8.2f}{m('errSR_pct'):8.2f}{m('errphi_pct'):8.2f}")

        print("\nSaturation pressure (bar) from phi_L = phi_V")
        print(f"{'gas':15s}{'T':>8s}{'ref':>9s}" + "".join(f"{e:>9s}" for e in EOS_PARAMS))
        for r in psat_rows:
            print(f"{r['gas']:15s}{r['T_K']:8.2f}{r['Psat_ref_bar']:9.3f}"
                  + "".join(f"{r[f'Psat_{e}_bar']:9.3f}" for e in EOS_PARAMS))
    return rows, psat_rows


if __name__ == "__main__":
    run()
