# Cubic Equation of State Calculator — PVT, H, S and Fugacity of Real Gases

Thermodynamics computer project (Prof. S.K. Gupta), IIT Delhi, 2026.

**Submitted by:** Krishna Kumar (2025CH72107) · Tutorial Group 02

## What it does

For a pure gas at a given temperature and pressure, using the **van der Waals, Redlich–Kwong, Soave–Redlich–Kwong or Peng–Robinson** equation of state, the program computes:

| Part | Output |
|---|---|
| PVT | compressibility factor Z (all real roots found, stable root chosen), molar volume V |
| HW: H and S | residual (departure) enthalpy H<sup>R</sup> and entropy S<sup>R</sup>, and total H and S (reference: ideal gas at 298.15 K, 1 bar) |
| Project: fugacity | fugacity coefficient φ and fugacity f = φP |
| Bonus | saturation pressure from the equal-fugacity condition φ<sup>L</sup> = φ<sup>V</sup> |

Built-in gases: methane, ethane, propane, n-butane, nitrogen, oxygen, carbon dioxide, ammonia, water.

## Live website

**https://krishnakr9861.github.io/CHL2002-PROJECT/**: the same calculator running in the browser (`index.html`, a JavaScript port of `eos.py` that gives identical results). Enter gas, T, P and EOS to get Z, V, H, S, φ and f, compare all four EOS, find the saturation pressure and plot isotherms.

## How to run

The Python program is in the `python/` folder:

```bash
cd python
pip install -r requirements.txt

python main.py                         # interactive menu
python main.py methane 250 50 PR       # gas, T (K), P (bar), EOS
python validate.py                     # comparison with reference data -> results/
python make_plots.py                   # report figures -> figures/
```

Example output (`python main.py methane 250 50 PR`):

```
  Gas: methane    EOS: PR
  T = 250.00 K    P = 50.0000 bar    Tr = 1.3116    Pr = 1.0872
  PVT
    Z                = 0.81172
    Molar volume V   = 337.452 cm3/mol
  Thermodynamic properties  (ref: ideal gas, 298.15 K, 1 bar)
    H_R (departure)  = -1295.01 J/mol
    S_R (departure)  = -3.5908 J/(mol K)
    H                = -2907.55 J/mol
    S                = -42.0077 J/(mol K)
  Fugacity
    phi              = 0.82601
    f                = 41.3005 bar
```

## Files

| File | Contents |
|---|---|
| `index.html` | browser version of the calculator (the website) |
| `report.pdf` | 4-page project report |
| `python/data.py` | Tc, Pc, ω and ideal-gas Cp coefficients (Smith, Van Ness & Abbott, App. B & C) |
| `python/eos.py` | EOS parameters, cubic solver, H<sup>R</sup>, S<sup>R</sup>, ln φ, ideal-gas H and S, saturation pressure |
| `python/main.py` | command-line / menu interface |
| `python/validate.py` | compares Z, H<sup>R</sup>, S<sup>R</sup>, φ and P<sup>sat</sup> with CoolProp (NIST-quality reference EOS) |
| `python/make_plots.py` | figures used in the report |
| `python/requirements.txt` | Python libraries needed |

`main.py` needs only NumPy. CoolProp and Matplotlib are needed only for validation and plots.
