"""
Cubic EOS calculator - PVT, enthalpy, entropy and fugacity of pure real gases.

Run:
    python main.py                      (interactive menu)
    python main.py methane 250 50 PR    (gas, T in K, P in bar, EOS)
"""

import sys

from data import COMPONENTS
from eos import EOS_PARAMS, state, saturation_pressure


def print_state(s):
    print()
    print(f"  Gas: {s['gas']}    EOS: {s['eos']}")
    print(f"  T = {s['T']:.2f} K    P = {s['P'] / 1e5:.4f} bar    "
          f"Tr = {s['Tr']:.4f}    Pr = {s['Pr']:.4f}")
    print("  " + "-" * 60)
    print(f"  EOS constants      a = {s['a']:.5e} Pa m6/mol2,  b = {s['b']:.5e} m3/mol")
    print(f"  Dimensionless      A = {s['A']:.5f},  B = {s['B']:.5f},  real roots = {s['n_roots']}")
    print("  PVT")
    print(f"    Z                = {s['Z']:.5f}")
    print(f"    Molar volume V   = {s['V'] * 1e6:.3f} cm3/mol")
    print("  Thermodynamic properties  (ref: ideal gas, 298.15 K, 1 bar)")
    print(f"    H_R (departure)  = {s['HR']:.2f} J/mol")
    print(f"    S_R (departure)  = {s['SR']:.4f} J/(mol K)")
    print(f"    H                = {s['H']:.2f} J/mol")
    print(f"    S                = {s['S']:.4f} J/(mol K)")
    print("  Fugacity")
    print(f"    phi              = {s['phi']:.5f}")
    print(f"    f                = {s['f'] / 1e5:.4f} bar")
    print()


def ask(prompt, cast=str, default=None):
    txt = input(f"{prompt}" + (f" [{default}]" if default is not None else "") + ": ").strip()
    if not txt and default is not None:
        return cast(default)
    return cast(txt)


def interactive():
    print("=" * 64)
    print("  CUBIC EQUATION OF STATE CALCULATOR")
    print("  PVT  |  H and S  |  pure-component fugacity")
    print("=" * 64)
    gases = list(COMPONENTS)
    eoses = list(EOS_PARAMS)
    while True:
        print("\nGases:", ", ".join(f"{i + 1}.{g}" for i, g in enumerate(gases)))
        print("1) Properties at T, P     2) Saturation pressure     3) Compare all EOS     0) Quit")
        choice = ask("Choice", str, "1")
        if choice == "0":
            break
        try:
            g = gases[ask("Gas number", int, 1) - 1]
            T = ask("Temperature (K)", float, 300)
            if choice == "2":
                eos = ask(f"EOS {eoses}", str, "PR")
                Psat, ZL, ZV = saturation_pressure(g, T, eos)
                print(f"\n  Psat = {Psat / 1e5:.4f} bar   Z_liq = {ZL:.5f}   Z_vap = {ZV:.5f}\n")
                continue
            P = ask("Pressure (bar)", float, 10) * 1e5
            if choice == "3":
                print(f"\n  {'EOS':5s}{'Z':>10s}{'V cm3/mol':>12s}{'H_R J/mol':>12s}"
                      f"{'S_R J/molK':>12s}{'phi':>9s}")
                for e in eoses:
                    s = state(g, T, P, e)
                    print(f"  {e:5s}{s['Z']:10.5f}{s['V'] * 1e6:12.2f}{s['HR']:12.1f}"
                          f"{s['SR']:12.3f}{s['phi']:9.4f}")
                print()
            else:
                eos = ask(f"EOS {eoses}", str, "PR")
                phase = ask("Phase (auto/vapor/liquid)", str, "auto")
                print_state(state(g, T, P, eos, phase))
        except Exception as exc:  # keep the menu alive on bad input
            print(f"  Error: {exc}")


if __name__ == "__main__":
    if len(sys.argv) >= 4:
        eos = sys.argv[4] if len(sys.argv) > 4 else "PR"
        print_state(state(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]) * 1e5, eos))
    else:
        interactive()
