"""
_style.py -- one visual convention for all nine figures.

The figures were written at different points in the project and drifted apart:
different fonts, line widths, and three different names for the same parameter.
Importing this module at the top of a figure script fixes typography and line
weights; the dictionaries below fix colour and wording, and have to be used
explicitly.

    from _style import apply, C, PHASE, CLASS, LABEL

    apply()                       # rcParams
    ax.plot(x, y, color=C["xor"])
    ax.set_xlabel(LABEL["w_hl"])

Conventions that matter for consistency with the text:

* Developmental phase is always builder -> switch -> governor, in that order,
  with the colours of PHASE, and the builder always on the left or at the top.
* Outcome classes use the wording of CLASS, which matches the classification
  defined in analysis/rescore_runs.py and quoted in the preprint: "failure to
  teach", "asymmetric", "silent", "XOR", "OR".  Not "collapse", which was the
  label of the superseded absolute-rate classifier and survives only where the
  text explicitly discusses that older scheme.
* The hub-to-lateral weight is always "hub-to-lateral coupling", never "shunt
  strength" or "hub strength": it is one conductance that acts as depolarizing
  drive in the builder phase and as inhibition in the governor phase, so naming
  it after either role alone is misleading.
"""
import matplotlib as mpl

# ---------------------------------------------------------------- palette
C = {
    # developmental phases
    "builder":   "#e08a00",
    "switch":    "#9e9e9e",
    "governor":  "#2e7d32",
    # outcome classes
    "xor":        "#2e7d32",
    "or":         "#e08a00",
    "asymmetric": "#7b5ea7",
    "not_taught": "#9e9e9e",
    "silent":     "#616161",
    # circuit elements (schematics)
    "input":      "#555555",
    "lateral":    "#1f4e79",
    "excitatory": "#2e7d32",
    "inhibitory": "#a11a1a",
    "depolar":    "#e08a00",
    "readout":    "#5a3a8a",
    # quantities
    "weight":     "#8a2be2",
    "rate":       "#c0392b",
    "ecl":        "#a11a1a",
}

PHASE = ("builder", "switch", "governor")

CLASS = {
    "not_taught": "failure to teach",
    "asymmetric": "asymmetric",
    "silent":     "silent",
    "XOR":        "XOR",
    "OR":         "OR",
}

LABEL = {
    "w_hl":      "hub$\\to$lateral coupling  $w_{\\mathrm{HL}}$ (nS per synapse)",
    "t_build":   "builder duration  $t_{\\mathrm{build}}$ (s)",
    "ecl":       "$E_{Cl}$ (mV)",
    "arm":       "arm weight  $w$",
    "xor_index": "xor\\_index",
    "rate":      "read-out rate (Hz)",
    "r11":       "coincidence rate  $R_{11}$ (Hz)",
    "cliff":     "firing cliff",
}

# Values quoted in the text, so figures and prose cannot drift apart.
CLIFF = 0.21
XI_XOR = 0.5
VT_MV = -50.0
EL_MV = -65.0


def apply(scale=1.0):
    """Set rcParams. Call once, before creating any figure."""
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10 * scale,
        "axes.titlesize": 11.5 * scale,
        "axes.labelsize": 11 * scale,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.9,
        "xtick.labelsize": 9.5 * scale,
        "ytick.labelsize": 9.5 * scale,
        "legend.fontsize": 9.5 * scale,
        "legend.frameon": False,
        "lines.linewidth": 2.2,
        "lines.markersize": 6,
        "figure.dpi": 160,
        "savefig.dpi": 160,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
    })


def phase_bands(ax, t_build, t_ramp, t_gov, alpha=0.55):
    """Shade the three developmental phases on a time axis, always in order."""
    ax.axvspan(0, t_build, color=C["builder"], alpha=alpha * 0.35, zorder=0)
    ax.axvspan(t_build, t_build + t_ramp, color=C["switch"],
               alpha=alpha * 0.45, zorder=0)
    ax.axvspan(t_build + t_ramp, t_build + t_ramp + t_gov,
               color=C["governor"], alpha=alpha * 0.30, zorder=0)
