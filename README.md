# EmergingXOR

Code, data and figures for the preprint

> **From homeostasis to computation: could the GABA switch build the core of the signed-XOR motif?**
> Jesús Marco de Lucas (IFCA, CSIC–Universidad de Cantabria).

The paper asks how a developmental **GABA switch** — the maturation of chloride
extrusion that flips GABA from depolarizing to shunting — could turn a
pre-positioned inhibitory hub from a **builder** (depolarizing coactivator that
lifts subthreshold excitatory arms over the firing cliff) into a **governor**
(shunting, coincidence-suppressing gain control), thereby maturing the
coincidence-suppression core of a signed-XOR motif. All quantitative claims are
backed by conductance-based [Brian2](https://briansimulator.org) simulations.

## Headline result

Over **100 random seeds** the full builder→switch→governor sequence produces an
XOR-like read-out **deterministically** (XOR 100/100, both arms grown 100/100,
median `xor_index` 0.78). A single deleted phase abolishes it (collapse / OR),
and a fixed chloride transition with plasticity frozen (`eta=0`) builds nothing.
The one thing coincidence suppression depends on is **hub input correlation**:
as hub afferents become independent (`shared_frac` 1.0 → 0.5) the read-out
degrades smoothly from XOR to OR, while the arms still always grow. The GABA
switch therefore **matures a preconfigured XOR scaffold** rather than carving
XOR de novo. See `robustness_sweep_en.png`.

## Repository layout

```
EmergingXOR.tex / .pdf        the preprint (figures referenced by bare filename)
references.bib                bibliography
app_growth_robustness.tex     drop-in robustness paragraph + caption for the appendix

phase_sim3.py                 the two-phase conductance-based model (canonical)
window2d.py / window_sweep.py 2D / 1D developmental-window sweeps
run_phase_experiments.py      runs phase_sim3 -> phase_factorial.json + phase_traces2.npz
run_robustness.py             multi-seed robustness battery (parallel) -> data/robust_*.json

make_motif_mature.py          Fig: mature signed-XOR motif        (_helpers.py)
make_motif_dev.py             Fig: developmental sequence         (_helpers.py)
make_preprint_fig.py          Fig: builder->governor curve
demo_brian2.py                Fig: the GABA switch (runs its own sim)
make_phase_figures.py         Figs: trace, factorial, window (two-phase)
make_robustness_fig.py        Fig: hub-independence robustness sweep

*_en.png                      the 8 figures used by the paper
data/*.json                   sweep results (window2d + robustness) and the
                              hardened-column spot-check
```

## Requirements

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```
(`brian2`, `numpy`, `matplotlib`, `pillow`.) The simulations use Brian2's
`numpy` codegen target for portability/determinism (works identically on x86-64
and ARM64/Grace); switch to `cython` at the top of `phase_sim3.py` for speed
where your build supports it.

## Reproducing the figures

```bash
# static schematics + the switch simulation
python3 make_motif_mature.py      # motif_mature_en.png
python3 make_motif_dev.py         # motif_dev_en.png
python3 make_preprint_fig.py      # preprint_fig_en.png
python3 demo_brian2.py            # demo_brian2_en.png

# developmental-window sweep (slow; a copy of the result ships in data/)
python3 window2d.py --seeds 3 --workers 8      # -> window2d.json

# two-phase data + figures (a few minutes)
python3 run_phase_experiments.py               # -> phase_factorial.json, phase_traces2.npz
python3 make_phase_figures.py                  # phase_trace_en / phase_factorial_en / window2d_en

# robustness battery + figure (results ship in data/)
for sf in 1.0 0.95 0.9 0.85 0.8 0.7 0.6 0.5; do
  python3 run_robustness.py --seeds 100 200 --shared $sf --workers 16 --t_build 18 --t_gov 12
done
python3 make_robustness_fig.py                 # robustness_sweep_en.png
```

`run_robustness.py` classifies **each** run individually and reports the XOR
fraction, both-grown fraction, and the median/IQR of the read-out rates and of
the continuous indices `xor_index = (min(R10,R01) − R11)/min(R10,R01)` and
`symmetry`.

## License

Code, data and figures are released under **CC BY-NC 4.0** (Creative Commons
Attribution–NonCommercial 4.0 International). See `LICENSE`.

## Citation

If you use this material, please cite the preprint (see `EmergingXOR.pdf`).
