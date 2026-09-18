"""
Supplementary code for:
"A Minimal Transport Model Reconciles Conflicting Reports of Oatp1a4 Zonation
in Rat Liver" (CPT:PSP Perspective submission)

A 24-compartment, single-pass sinusoidal transport model of the rat hepatic
lobule (portal vein -> central vein), used to compare three candidate spatial
profiles for Oatp1a4 activity against the SR-101 imaging data and kinetic
parameters reported in Akanuma et al., Drug Metab Pharmacokinet 2019;34:239-246.

No external dependencies (pure Python). Reproduces Figure 1A and 1B.
"""

N = 24                      # number of well-mixed compartments, PV -> CV
xs = [(i + 0.5) / N for i in range(N)]   # lobule position, 0 = PV, 1 = CV

# Kinetic parameters
KM = 8.2      # uM, high-affinity Km measured in isolated rat hepatocytes (Akanuma et al. 2019)
KI = 0.3      # uM, illustrative digoxin inhibition constant
F = 8.0       # sinusoidal blood flow (arbitrary units / min)
VB = 0.15     # fractional blood volume per compartment
VH = 0.6      # fractional hepatocyte volume per compartment
KEF = 0.8     # first-order basolateral efflux rate (/min)
VMAX0 = 10.0  # total hepatic Oatp1a4 capacity (arbitrary units), held equal across hypotheses
DT = 0.0025   # integration time step (min)


def vmax_profile(hypothesis):
    """Spatial Vmax(x) profile for each zonation hypothesis, normalized so
    all three hypotheses have identical total Oatp1a4 capacity."""
    if hypothesis == "H1":       # sharp pericentral (Akanuma et al. 2019, IHC)
        raw = [0.04 + 0.96 * x**3 for x in xs]
    elif hypothesis == "H2":     # non-zonated (van de Graaf et al. 2024)
        raw = [1.0 for _ in xs]
    elif hypothesis == "H3":     # continuous gradient (synthesis)
        raw = [0.2 + 0.8 * x for x in xs]
    else:
        raise ValueError(hypothesis)
    mean_raw = sum(raw) / len(raw)
    return [r / mean_raw * VMAX0 for r in raw]


def simulate(hypothesis, digoxin_uM=0.0, dose_uM=10.0, t_end_min=10.0):
    """Integrate the compartmental model forward in time (explicit Euler).

    Returns (Cb, Ch): blood and hepatocyte concentration arrays (length N,
    PV -> CV) at t = t_end_min.
    """
    vmax = vmax_profile(hypothesis)
    km_app = KM * (1 + digoxin_uM / KI) if digoxin_uM > 0 else KM

    Cb = [0.0] * N   # blood concentration per compartment
    Ch = [0.0] * N   # hepatocyte concentration per compartment
    nsteps = int(t_end_min / DT)

    for _ in range(nsteps):
        Cb_in = [dose_uM] + Cb[:-1]     # inflow to compartment i = outflow of i-1
        new_Cb, new_Ch = [0.0] * N, [0.0] * N
        for i in range(N):
            j_in = vmax[i] * Cb[i] / (km_app + Cb[i])   # saturable active influx
            j_out = KEF * Ch[i]                          # passive basolateral efflux
            dCb = F / VB * (Cb_in[i] - Cb[i]) - (j_in - j_out) / VB
            dCh = (j_in - j_out) / VH
            new_Cb[i] = max(0.0, Cb[i] + DT * dCb)
            new_Ch[i] = max(0.0, Ch[i] + DT * dCh)
        Cb, Ch = new_Cb, new_Ch

    return Cb, Ch


def cv_pv_ratio(Ch, n_edge=3):
    """Mean hepatocyte concentration over the outermost n_edge compartments
    at each end, expressed as a CV:PV ratio (Figure 1A summary statistic)."""
    pv = sum(Ch[:n_edge]) / n_edge
    cv = sum(Ch[-n_edge:]) / n_edge
    return cv / max(pv, 1e-9)


if __name__ == "__main__":
    print("Figure 1A -- baseline (10 uM SR-101, no digoxin, t = 10 min)")
    for h in ["H1", "H2", "H3"]:
        Cb, Ch = simulate(h, digoxin_uM=0.0, dose_uM=10.0, t_end_min=10.0)
        print(f"  {h}: CV:PV ratio = {cv_pv_ratio(Ch):.2f}")

    print()
    print("Figure 1B -- CV signal (% of no-digoxin baseline) vs. digoxin dose")
    digoxin_doses = [0, 0.3, 0.6, 1.0, 1.5, 2.0, 2.4, 3.0, 4.0, 5.0]
    for h in ["H1", "H2", "H3"]:
        cv0 = None
        row = []
        for d in digoxin_doses:
            _, Ch = simulate(h, digoxin_uM=d, dose_uM=10.0, t_end_min=10.0)
            cv = sum(Ch[-3:]) / 3
            cv0 = cv0 or cv
            row.append(f"{100 * cv / cv0:5.1f}%")
        print(f"  {h}: " + "  ".join(row))
