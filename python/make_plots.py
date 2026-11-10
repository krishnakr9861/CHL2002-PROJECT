import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from eos import EOS_PARAMS, state, saturation_pressure
from validate import reference

COLORS = {"vdW": "#2a78d6", "RK": "#eb6834", "SRK": "#1baf7a", "PR": "#e87ba4"}
REF_STYLE = dict(color="#0b0b0b", marker="o", s=18, zorder=5, label="Reference (NIST-quality EOS)")

plt.rcParams.update({
    "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#c3c2b7", "axes.grid": True, "grid.color": "#e1e0d9",
    "grid.linewidth": 0.6, "lines.linewidth": 1.6, "legend.frameon": False,
    "savefig.dpi": 200,
})


def property_panels():
    gas, T = "carbon dioxide", 350.0
    P = np.linspace(1, 200, 120)
    Pref = np.linspace(10, 200, 12)
    refs = [reference(gas, T, p * 1e5) for p in Pref]

    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.4))
    panels = [("Z", "Compressibility factor Z", 1),
              ("HR", "Residual enthalpy $H^R$ (kJ/mol)", 1e-3),
              ("SR", "Residual entropy $S^R$ (J/mol K)", 1),
              ("phi", r"Fugacity coefficient $\phi$", 1)]
    for ax, (key, label, scale) in zip(axs.flat, panels):
        for eos in EOS_PARAMS:
            vals = [state(gas, T, p * 1e5, eos)[key] * scale for p in P]
            ax.plot(P, vals, color=COLORS[eos], label=eos)
        ax.scatter(Pref, [r[key] * scale for r in refs], **REF_STYLE)
        ax.set_xlabel("P (bar)")
        ax.set_ylabel(label)
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, fontsize=8)
    fig.suptitle(f"CO$_2$ at {T:.0f} K: four cubic EOS vs reference data", fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig("figures/properties.png")
    plt.close(fig)


def isotherm_panels():
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9))

    # (a) Z vs Pr for methane, PR EOS, three isotherms
    ax = axs[0]
    temps = [200.0, 250.0, 350.0]
    tcolors = ["#2a78d6", "#eb6834", "#1baf7a"]
    P = np.linspace(1, 300, 150)
    for T, c in zip(temps, tcolors):
        ax.plot(P, [state("methane", T, p * 1e5, "PR")["Z"] for p in P], color=c, label=f"{T:.0f} K")
        Pr = np.linspace(20, 300, 8)
        ax.scatter(Pr, [reference("methane", T, p * 1e5)["Z"] for p in Pr],
                   color=c, edgecolor="#0b0b0b", linewidth=0.6, s=16, zorder=5)
    ax.set_xlabel("P (bar)")
    ax.set_ylabel("Z")
    ax.set_title("(a) Methane, PR (lines) vs reference (dots)", fontsize=9)
    ax.legend(fontsize=8)

    # (b) P-V isotherm of propane below Tc with the equal-fugacity Psat
    ax = axs[1]
    T = 340.0
    from data import R, get_component
    from eos import eos_constants
    comp = get_component("propane")
    a, b, _ = eos_constants(comp, T, "PR")
    V = np.geomspace(b * 1.15, 5e-3, 600)
    Pcurve = R * T / (V - b) - a / (V * V + 2 * b * V - b * b)
    ax.plot(V * 1e6, Pcurve / 1e5, color="#2a78d6", label="PR isotherm, 340 K")
    Psat, ZL, ZV = saturation_pressure("propane", T, "PR")
    VL, VV = ZL * R * T / Psat, ZV * R * T / Psat
    ax.hlines(Psat / 1e5, VL * 1e6, VV * 1e6, color="#eb6834", linestyle="--",
              label=f"$P^{{sat}}$ = {Psat / 1e5:.1f} bar ($f^L = f^V$)")
    ax.set_xscale("log")
    ax.set_ylim(-20, 80)
    ax.set_xlabel("V (cm$^3$/mol)")
    ax.set_ylabel("P (bar)")
    ax.set_title("(b) Propane: van der Waals loop & saturation", fontsize=9)
    ax.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig("figures/isotherms.png")
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    property_panels()
    isotherm_panels()
    print("Saved figures/properties.png and figures/isotherms.png")
