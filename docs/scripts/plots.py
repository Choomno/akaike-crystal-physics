"""Draw figures from verified canonical results; never duplicate fitted values."""
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from models import ROOT, debye, einstein

OUT = ROOT / "docs" / "img"
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


def save(fig, stem, theme):
    for extension in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{extension}", facecolor=theme["surface"],
                    dpi=180, metadata=({"Creator": "akaike-crystal-physics/docs/scripts/plots.py",
                                       "CreationDate": None, "ModDate": None} if extension == "pdf"
                                      else {"Software": "akaike-crystal-physics/docs/scripts/plots.py"}))
    plt.close(fig)


def main():
    audit = json.loads((ROOT / "docs/generated/audit.json").read_text())
    for relative, expected in audit["inputs"].items():
        if hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Input changed: {relative}; run verify.py first")
    OUT.mkdir(exist_ok=True)
    data = audit["datasets"]["swaminathan_a"]
    temperature, observed = np.array(data["x"]), np.array(data["y"])
    grid = np.linspace(temperature.min(), temperature.max(), 400)
    models = audit["urea_models"]
    labels = ["Linear", "Quadratic", "Cubic", "Einstein", "Debye"]
    for mode, theme in THEMES.items():
        fig, (ax, ranking) = plt.subplots(1, 2, figsize=(10.5, 4.6),
                                         gridspec_kw={"width_ratios": [1.25, 1]})
        fig.patch.set_facecolor(theme["surface"])
        style(ax, theme); style(ranking, theme)
        colors = [theme["series"][0], theme["series"][1], theme["muted"], theme["series"][2], "#ae70c9"]
        for i, record in enumerate(models["corrected"]):
            if i == 2:  # cubic is present in the complete candidate ranking at right
                continue
            beta = record["parameters"]
            fitted = (einstein(beta, grid) if i == 3 else debye(beta, grid) if i == 4
                      else np.polynomial.polynomial.polyval(grid, beta))
            ax.plot(grid, fitted, color=colors[i], lw=1.8,
                    ls=["--", ":", "-.", "-", "--"][i], label=labels[i])
        ax.plot(temperature, observed, "o", ms=5, mfc=theme["surface"], mec=theme["ink"],
                mew=1.2, label="Measured", zorder=5)
        ax.set(xlabel="Temperature (K)", ylabel="Lattice parameter a (Å)")
        ax.set_title("Independent refits", loc="left", color=theme["ink"], fontsize=11)
        legend = ax.legend(frameon=False, fontsize=8, loc="upper left")
        for text in legend.get_texts(): text.set_color(theme["ink2"])
        positions = np.arange(5)
        for offset, key, color, title in ((-.18, "archived", theme["series"][0], "2020 archive"),
                                         (.18, "corrected", theme["series"][1], "Independent refit")):
            values = [r["weight"] for r in models[key]]
            ranking.barh(positions + offset, values, height=.32, color=color, label=title)
            for position, value in zip(positions + offset, values):
                ranking.text(value + .012, position, f"{value:.3f}" if value >= .001 else "<0.001",
                             va="center", color=theme["ink2"], fontsize=8)
        ranking.set(yticks=positions, yticklabels=labels, xlabel="Akaike weight", xlim=(0, .76))
        ranking.invert_yaxis()
        ranking.set_title("Same five candidate models", loc="left", color=theme["ink"], fontsize=11)
        legend = ranking.legend(frameon=False, fontsize=8, loc="lower right")
        for text in legend.get_texts(): text.set_color(theme["ink2"])
        fig.suptitle("Urea a(T): correcting the Debye fit changes the ranking", x=.06, ha="left",
                     color=theme["ink"], fontsize=13)
        fig.text(.06, .025, "Swaminathan data: 7 observations, 12–173 K. Raw a fitted; K includes error variance.\n"
                 "Curves cover the observed interval. Weights compare this candidate set only; <0.001 is not zero.",
                 color=theme["ink2"], fontsize=8)
        fig.tight_layout(rect=(0, .12, 1, .94))
        save(fig, f"urea_a_{mode}", theme)

        anthracene = audit["datasets"]["anthracene_a"]
        degrees = np.arange(1, 11)
        fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True)
        fig.patch.set_facecolor(theme["surface"])
        for axis in axes: style(axis, theme)
        for key, color, offset, title in (("archived", theme["series"][0], -.17, "2020 archive"),
                                         ("corrected", theme["series"][1], .17, "Scaled least squares")):
            rows = anthracene[key]
            axes[0].bar(degrees+offset, [r["weight"] for r in rows], width=.31, color=color, label=title)
            axes[1].plot(degrees, [1-r["adjusted_r2"] for r in rows],
                         "o--" if key == "archived" else "s-", ms=4, lw=1.2, color=color, label=title)
        axes[0].set(ylabel="Akaike weight", ylim=(0, 1))
        legend=axes[0].legend(frameon=False, fontsize=9)
        for text in legend.get_texts(): text.set_color(theme["ink2"])
        axes[1].set(yscale="log", ylabel="1 − adjusted R² (log scale)", xlabel="Polynomial degree", xticks=degrees)
        axes[0].set_title("Anthracene a(P): archived and independently fitted models", loc="left",
                          color=theme["ink"], fontsize=11)
        fig.text(.12, .02, "18 observations, 1.05–22.7 GPa; all ten degrees retained in each ranking.\n"
                 "Historical normal equations lose accuracy at degrees 8–10. A numerically stable fit\n"
                 "does not establish that a high-degree polynomial is a suitable physical model.",
                 color=theme["ink2"], fontsize=8)
        fig.tight_layout(rect=(0, .14, 1, 1))
        save(fig, f"anthracene_a_{mode}", theme)
    print("Wrote four figures in PNG and PDF to docs/img/")


if __name__ == "__main__":
    main()
