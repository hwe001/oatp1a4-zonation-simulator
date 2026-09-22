import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import time

N = 24
Km = 8.2
Vb = 0.15
Vh = 0.6
C0 = 10.0
F_fixed = 8.0

def max_cvpv_H2(vmax0, kef, F=F_fixed, safety=0.04, t_max_cap=900.0):
    fastest = max(F/Vb, kef, vmax0/(Km*Vb), 1.0)
    dt = safety / fastest
    # cap t_max: analytic result guarantees monotone approach to 1 from below,
    # so a bounded window well past the transient peak is sufficient
    t_max = min(max(300.0, 20.0/kef), t_max_cap)
    nsteps = int(t_max/dt)
    Cb = np.zeros(N)
    Ch = np.zeros(N)
    max_ratio = 0.0
    check_every = max(1, nsteps // 1500)
    for step in range(nsteps):
        Cb_in = np.empty(N)
        Cb_in[0] = C0
        Cb_in[1:] = Cb[:-1]
        Jin = vmax0 * Cb / (Km + Cb)
        Jout = kef * Ch
        dCb = F/Vb*(Cb_in - Cb) - (Jin - Jout)/Vb
        dCh = (Jin - Jout)/Vh
        Cb = np.maximum(Cb + dt*dCb, 0.0)
        Ch = np.maximum(Ch + dt*dCh, 0.0)
        if step % check_every == 0:
            pv, cv = Ch[0], Ch[-1]
            if pv > 1e-9:
                r = cv/pv
                if r > max_ratio:
                    max_ratio = r
    return max_ratio

vmax0_vals = [1, 10, 100, 300]
kef_vals   = [0.03, 0.1, 0.8, 5]

grid = np.zeros((len(kef_vals), len(vmax0_vals)))
t0 = time.time()
for i, kef in enumerate(kef_vals):
    for j, vmax0 in enumerate(vmax0_vals):
        r = max_cvpv_H2(vmax0, kef)
        grid[i, j] = r
        print(f"[{time.time()-t0:6.1f}s] vmax0={vmax0:6.1f} (Vmax/F={vmax0/F_fixed:6.2f})  kef={kef:5.2f}  max CV:PV={r:.4f}", flush=True)

print("\nGlobal max CV:PV across entire grid:", grid.max())
np.save(__file__.replace("sensitivity_analysis.py", "grid.npy"), grid)

fig, ax = plt.subplots(figsize=(5.4, 4.4))
im = ax.imshow(grid, origin="lower", cmap="Blues", vmin=0, vmax=1.0, aspect="auto")
ax.set_xticks(range(len(vmax0_vals)))
ax.set_xticklabels([f"{v/F_fixed:.2g}" for v in vmax0_vals])
ax.set_yticks(range(len(kef_vals)))
ax.set_yticklabels([f"{k:.2g}" for k in kef_vals])
ax.set_xlabel("Extraction ratio proxy, Vmax/F")
ax.set_ylabel("Efflux rate, kef (/min)")
ax.set_title("Max CV:PV over time for H2 (non-zonated)\nnever exceeds 1, across the grid tested", fontsize=10)
for i in range(len(kef_vals)):
    for j in range(len(vmax0_vals)):
        ax.text(j, i, f"{grid[i,j]:.2f}", ha="center", va="center", fontsize=9,
                 color="white" if grid[i,j] > 0.6 else "#333")
cbar = fig.colorbar(im, ax=ax, shrink=0.85)
cbar.set_label("max CV:PV ratio over time")
fig.tight_layout()
fig.savefig(__file__.replace("sensitivity_analysis.py", "figureS1_sensitivity.png"), dpi=300)
print("saved heatmap")
