"""The figures of the paper, drawn from the data that ../figure_data.jl writes to this folder (domain.json,
total_loss.json).

    python plot_figures.py

    domain.pdf        the domain D of (24) at theta = 0.1129, its Ranges 1-6 (Table 1), the discarded regions
                      D1, D2, D4, D6, D8 and the region R (26)
    total_loss.pdf    the total loss of each decomposition of Table 3 for 0.112 <= theta <= 0.118
"""
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True, "axes.unicode_minus": False, "font.size": 10,
    "pdf.fonttype": 42, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})
GRAY_D, GRAY_DISCARD, BLUE_R, INK = "#e8e7e3", "#9c9a92", "#2a78d6", "#52514e"
LINE_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]

# Positions of the labels, in (alpha, beta); a third entry is the point an arrow from the label points to.
RANGE_LABELS = {1: (0.443, 0.300), 2: (0.360, 0.337), 3: (0.312, 0.330), 4: (0.282, 0.305), 5: (0.235, 0.255),
                6: (0.128, 0.098, (0.140, 0.122))}
REGION_LABELS = {"D1": (0.475, 0.190), "D2": (0.362, 0.245), "D4": (0.312, 0.245), "D6": (0.300, 0.105, (0.284, 0.215)),
                 "D8": (0.256, 0.234), "R": (0.405, 0.112, (0.375, 0.136))}


def fill(ax, region, color):
    ax.fill_between(region["alpha"], region["lower"], region["upper"], facecolor=color, edgecolor=color, linewidth=0.3)


def outline(ax, pieces):
    """The outer boundary of the union of the pieces, which are ordered in alpha and adjacent."""
    a = [x for p in pieces for x in p["alpha"]]
    top = [x for p in pieces for x in p["upper"]]
    bottom = [x for p in pieces for x in p["lower"]]
    ax.plot(a + a[::-1] + a[:1], top + bottom[::-1] + top[:1], color="black", linewidth=0.6)


def label(ax, text, spec, **kw):
    if len(spec) == 2:
        ax.text(*spec, text, ha="center", va="center", **kw)
    else:
        ax.annotate(text, xy=spec[2], xytext=spec[:2], ha="center", va="center",
                    arrowprops=dict(arrowstyle="-", color="black", linewidth=0.5, shrinkA=1, shrinkB=0), **kw)


def domain_figure():
    data = json.loads((HERE / "domain.json").read_text())
    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    for piece in data["D"]:
        fill(ax, piece, GRAY_D)
    for region in data["discarded"]:
        fill(ax, region, GRAY_DISCARD)
    fill(ax, data["R"], BLUE_R)
    outline(ax, data["D"])
    for a, b0, b1 in data["range_verticals"]:
        ax.plot([a, a], [b0, b1], color=INK, linestyle=(0, (4, 3)), linewidth=0.7)
    ax.plot(data["range_diagonal"]["alpha"], data["range_diagonal"]["beta"], color=INK, linestyle=(0, (4, 3)),
            linewidth=0.7)
    for k, spec in RANGE_LABELS.items():
        label(ax, str(k), spec, color=INK, fontsize=11)
    for name, spec in REGION_LABELS.items():
        text = r"$\mathsf{R}$" if name == "R" else rf"$\mathsf{{D}}_{{{name[1:]}}}$"
        label(ax, text, spec, fontsize=11)
    ax.set_xlim(0.11, 0.51)
    ax.set_ylim(0.08, 0.35)
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$\beta$")
    fig.tight_layout()
    return fig


def total_loss_figure():
    data = json.loads((HERE / "total_loss.json").read_text())
    theta = data["theta"]
    fig, ax = plt.subplots(figsize=(5.4, 3.3))
    ax.axhline(1.0, color=INK, linestyle=(0, (4, 3)), linewidth=0.7)
    for color, curve in zip(LINE_COLORS, data["curves"]):
        ax.plot(theta, curve["total"], color=color, linewidth=1.5,
                label=re.sub(r"\b([RGPL])\b", r"$\\mathsf{\1}$", curve["label"]))
        ax.plot([curve["root"]], [1.0], "o", color=color, markersize=5, markeredgecolor="white", markeredgewidth=1)
    ax.set_xlim(theta[0], theta[-1])
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel("total loss")
    ax.grid(axis="y", color=GRAY_D, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    import sys

    preview = Path(sys.argv[1]) if len(sys.argv) > 1 else None      # optional directory for PNG previews
    for name, make in (("domain.pdf", domain_figure), ("total_loss.pdf", total_loss_figure)):
        fig = make()
        fig.savefig(HERE / name)
        if preview is not None:
            fig.savefig(preview / name.replace(".pdf", ".png"), dpi=200)
        plt.close(fig)
        print("wrote", HERE / name)
