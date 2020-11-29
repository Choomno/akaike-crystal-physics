# README figures, regenerated in 2026 from the 2020 fitted values / spreadsheets.
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
OUT = os.path.join(R, "docs", "img")

THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781",
                  grid="#e1e0d9", axis="#c3c2b7", series=["#2a78d6", "#eb6834", "#1baf7a"]),
    "dark":  dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
                  grid="#2c2c2a", axis="#383835", series=["#3987e5", "#d95926", "#199e70"]),
}

def style(ax, t):
    ax.set_facecolor(t["surface"])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(t["axis"])
    ax.tick_params(colors=t["muted"], labelcolor=t["ink2"], labelsize=9)
    ax.grid(True, color=t["grid"], linewidth=0.6)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(t["ink2"]); ax.yaxis.label.set_color(t["ink2"])

# ---------- Figure 1: urea lattice parameter a (Swaminathan data) ----------
sw = np.loadtxt(os.path.join(R, "data", "swaminathan.txt"))
T, a = sw[0], sw[1]
Tf = np.linspace(1, 185, 400)
b1 = np.polyfit(T, a, 1); b2 = np.polyfit(T, a, 2)
ein = 5.56491 + 0.0373932 / (np.exp(133.132 / Tf) - 1)          # 2020 nonlinear.xls
curves = [  # (label, y, Akaike weight from 2020 AIC.xls, dash)
    ("Linear", np.polyval(b1, Tf), 0.125, (0, (5, 3))),
    ("Quadratic", np.polyval(b2, Tf), 0.254, (0, (1.5, 2))),
    ("Einstein (3 params)", ein, 0.611, "solid"),
]

for mode, t in THEMES.items():
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    fig.patch.set_facecolor(t["surface"]); style(ax, t)
    for (lab, y, w, ls), col, dy in zip(curves, t["series"], (-8, 7, -2)):
        ax.plot(Tf, y, color=col, lw=2, ls=ls, label=f"{lab}  ·  w = {w:.2f}")
        ax.annotate(lab, (Tf[-1], y[-1]), xytext=(6, dy), textcoords="offset points",
                    va="center", fontsize=9, color=t["ink2"])
    ax.plot(T, a, "o", ms=7, mfc=t["surface"], mec=t["ink"], mew=1.6, label="Measured (Swaminathan et al.)", zorder=5)
    ax.set_xlim(0, 215)
    ax.set_xlabel("Temperature (K)"); ax.set_ylabel("Lattice parameter a (Å)")
    ax.set_title("Urea, lattice parameter a: candidate models and Akaike weights",
                 loc="left", fontsize=11, color=t["ink"])
    leg = ax.legend(loc="upper left", frameon=False, fontsize=9)
    for txt in leg.get_texts(): txt.set_color(t["ink2"])
    fig.text(0.01, 0.01, "w = Akaike weight from the 2020 run (cubic ≈ 0; Debye 0.01, see Known issues). "
             "Figure regenerated 2026 from the 2020 fitted values.", fontsize=7, color=t["muted"])
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(OUT, f"urea_a_{mode}.png"), facecolor=t["surface"])
    plt.close(fig)

# ---------- Figure 2: anthracene a, Akaike weight vs. adjusted R^2 ----------
deg = np.arange(1, 10)
w = np.array([0, 0, 3.3788e-06, 0.0106691, 0.201719, 0.734541, 0.0528413, 0.000226343, 1.40049e-07])
r2 = np.array([0.905174, 0.986378, 0.99728, 0.999075, 0.999469, 0.999656, 0.999684, 0.99965, 0.999601])

for mode, t in THEMES.items():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.4), dpi=150, sharex=True,
                                   gridspec_kw=dict(hspace=0.35))
    fig.patch.set_facecolor(t["surface"])
    col = t["series"][0]
    style(ax1, t); style(ax2, t)
    ax1.bar(deg, w, width=0.6, color=col)
    ax1.set_ylabel("Akaike weight"); ax1.set_ylim(0, 1)
    ax1.set_title("Akaike weight picks one model (degree 6, w = 0.73)",
                  loc="left", fontsize=10, color=t["ink"])
    ax1.annotate("0.73", (6, w[5]), xytext=(0, 4), textcoords="offset points",
                 ha="center", fontsize=9, color=t["ink2"])
    ax2.plot(deg, 1 - r2, "-o", color=col, lw=2, ms=6, mfc=col, mec=t["surface"], mew=2)
    ax2.set_yscale("log")
    ax2.set_ylabel("1 − adjusted R²\n(log scale)")
    ax2.set_title("Adjusted R² creeps up to degree 7 (0.99908 → 0.99968), with no clear stopping point",
                  loc="left", fontsize=10, color=t["ink"])
    ax2.set_xlabel("Polynomial degree"); ax2.set_xticks(deg)
    fig.suptitle("Anthracene under pressure, lattice parameter a (Oehzelt et al. data)",
                 x=0.01, ha="left", fontsize=11, color=t["ink"])
    fig.text(0.01, 0.01, "Degree 10 omitted (numerically unstable fit, see Known issues). "
             "Values from the 2020 output; figure regenerated 2026.", fontsize=7, color=t["muted"])
    fig.subplots_adjust(left=0.14, right=0.97, top=0.88, bottom=0.12)
    fig.savefig(os.path.join(OUT, f"anthracene_a_{mode}.png"), facecolor=t["surface"])
    plt.close(fig)
print("done")
