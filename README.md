# EmergingXOR

Code and data for **"From homeostasis to computation: how the GABA switch could
mature the core of a signed-XOR motif"** (Marco de Lucas, preprint, 2026).

Every figure in the preprint can be regenerated from this repository. Nothing is
plotted from intermediate files of unclear provenance: each figure script reads
either nothing at all or a cache written by a named simulation script in
`simulations/`.

---

## Quick start

```bash
conda env create -f environment.yml       # or: pip install -r requirements.txt
conda activate emergingxor

python reproduce_all.py --quick           # schematics + fast figures, ~10 min
python reproduce_all.py                   # everything, several hours
```

`--quick` skips the three long sweeps (Figs. 4, 7, 8) and uses whatever caches
are already in `data/`. Run it first: if it fails, the environment is wrong and
there is no point launching the long jobs.

---

## What produces what

| Figure | Script | Input | Runtime | Output |
|---|---|---|---|---|
| Fig. 1 | `figures/tikz/fig1_mature_motif.tex` | none (TikZ, inlined in the manuscript) | — | vector, no PNG |
| Fig. 2 | `figures/fig2_development.py` | none (schematic) | seconds | `motif_dev_en.png` |
| Fig. 3 | `simulations/gaba_switch.py` | none | ~2 min | `demo_brian2_en.png` |
| Fig. 4 | `figures/fig4_ecl_setpoints.py` | `data/fig4/ecl_sweep.json` | seconds | `ecl_sweep_en.png` |
| ↳ its data | `simulations/ecl_sweep.py --seeds 20` | none | ~4 h (8 cores) | `data/fig4/ecl_sweep.json` |
| Figs. 5–6 | `figures/fig5_fig6_growth_controls.py` | `data/fig5_6/` | seconds | `phase_trace_en.png`, `phase_factorial_en.png` |
| ↳ its data | `simulations/phase_controls.py` | none | ~15 min | `phase_factorial.json`, `phase_traces2.npz` |
| Fig. 7 | `simulations/developmental_window.py --plot` | `data/fig7/window2d.json` | seconds | `window2d_en.png` |
| ↳ its data | `simulations/developmental_window.py --seeds 20` | none | ~8 h (8 cores) | `data/fig7/window2d.json` |
| Fig. 8 | `figures/fig8_hub_correlation.py` | `data/fig8/robust_*.json` | seconds | `robustness_sweep_en.png` |
| ↳ its data | `simulations/hub_correlation.py` (sweep, see below) | none | ~6 h (8 cores) | `data/fig8/robust_shared*.json` |
| Fig. 9 | `figures/fig9_builder_governor.py` | none (schematic) | seconds | `preprint_fig_en.png` |

Figure 1 is drawn in TikZ and inlined in the manuscript, so it is vector rather
than raster and has no generating script in `figures/`. To preview or edit it
without recompiling the whole paper:

```bash
cd figures/tikz && pdflatex standalone.tex
```

`standalone.tex` `\input`s the same file the manuscript includes, so the preview
cannot drift from the published figure.

The hub-correlation sweep of Fig. 8 is a loop over the shared-afferent fraction:

```bash
for sf in 1.0 0.95 0.9 0.8 0.7 0.6 0.5; do
  python simulations/hub_correlation.py --seeds 100 140 --shared $sf \
      --workers 8 --t_build 18 --t_gov 12 --out data/fig8/robust_shared${sf}.json
done
```

Runtimes assume 8 worker processes with Brian2's `numpy` code-generation target.
Do not switch to `cython`: the workers compile concurrently into a shared cache
and collide.

---

## Numbers quoted in the text

| Claim | Where | How to check |
|---|---|---|
| `g_inh(2)/g_inh(1) = 5.9`, excess 171 nS vs σ_eff ≈ 17 nS | §A.3, Prediction 1 | `python analysis/measure_supralinearity.py` |
| Builder-only arms compute XOR at mature E_Cl (`xor_index = 0.62`) | §A.3, Fig. 6 caption | same script, over-growth block |
| Re-scoring moves 161/800 labels | §A.3 | `python analysis/rescore_runs.py --file data/fig7/window2d.json` |
| Reference DOIs | bibliography | `python analysis/check_dois.py --file ../EmergingXOR.tex` |

---

## Classification

Every classified figure uses the same rule, defined once in
`analysis/rescore_runs.py` and duplicated (deliberately, so a sweep can label
runs as it goes) in `simulations/developmental_window.py::label`. The two are
kept in sync; `reproduce_all.py --check` asserts they agree on the cached runs.

With `w_I, w_P` the two arm weights, cliff `= 0.21`, and

```
xor_index = (min(R10, R01) − R11) / min(R10, R01)
```

a run is labelled:

- **failure to teach** — neither arm cleared the cliff, `max(w_I, w_P) ≤ 0.21`
- **asymmetric** — exactly one arm cleared it
- **silent** — both cleared it but `min(R10, R01) = 0`
- **XOR** — `xor_index ≥ 0.5`
- **OR** — otherwise

Read-out rates in hertz are reported as descriptors and never enter the
classification. An earlier version of this code gated on an absolute rate
(`min(R10,R01) < 25 Hz → collapse`), which mislabelled low-rate but selective
circuits; §A.3 of the preprint explains the change, and `rescore_runs.py`
reproduces both labellings so they can be compared.

---

## Renamed from the working version

If you are following the preprint's development history, the scripts were
renamed on release:

| was | is now |
|---|---|
| `demo_brian2.py` | `simulations/gaba_switch.py` |
| `phase_sim3.py` | `simulations/two_phase_model.py` |
| `window2d.py` | `simulations/developmental_window.py` |
| `run_robustness.py` | `simulations/hub_correlation.py` |
| `run_phase_experiments.py` | `simulations/phase_controls.py` |
| `rescore.py` | `analysis/rescore_runs.py` |
| `make_motif_dev.py` | `figures/fig2_development.py` |
| `make_motif_mature.py` | superseded by `figures/tikz/fig1_mature_motif.tex` |
| `make_ecl_fig.py` | `figures/fig4_ecl_setpoints.py` |
| `make_phase_figures.py` | `figures/fig5_fig6_growth_controls.py` |
| `make_robustness_fig.py` | `figures/fig8_hub_correlation.py` |
| `make_preprint_fig.py` | `figures/fig9_builder_governor.py` |

Two scripts from the working version are **not** included: `demo_grow_xor.py`
and `demo_grow_xor_bcm.py`. They are an earlier single-phase model whose arms
start above the firing cliff, and their behaviour contradicts the two-phase
result reported in the preprint. They were superseded, not used for any figure.

---

## Parameters

The model parameters live in `simulations/two_phase_model.py` and match the
preprint exactly: builder `E_Cl = −40 mV`, ramp to `−65 mV` over 4 s, governor
`E_Cl = −65 mV`; `w₀ = 0.10`, cliff ≈ 0.21, `w ∈ [0.05, 2.6]`; `τ_pre = 15 ms`,
`τ_Ca = 60 ms`, `κ = 9`, `τ_θ = 800 ms`, `η = 3×10⁻⁵`; 20 hub cells;
`V_T = −50 mV`, `E_L = −65 mV`. The network of Fig. 3
(`simulations/gaba_switch.py`) uses `−45 mV` and `−75 mV`, as stated in §A.2.

Seeds are generated, not hard-coded: every sweep takes `--seed0` (default 11)
and `--seeds N`, and uses `range(seed0, seed0 + N)`. The seed of each run is
stored in its cache entry, so any single run can be repeated in isolation.

---

## Licence

Everything here -- manuscript, code, figures and data -- is released under
**CC BY-NC 4.0** (Attribution-NonCommercial). See `LICENSE`.

Two things that licence does *not* do, and which matter here:

* It covers copyright only. A European patent application, *Signed XOR
  Feedback* (No. 26382252.0), has been filed on the mechanism this code
  implements, and no patent rights are granted by the licence. If you intend to
  practise the mechanism, for any purpose, contact the author.
* "NonCommercial" is not defined precisely in the licence text. If your
  intended use is anywhere near that boundary -- industry affiliation,
  industrially funded research, a paid service -- ask rather than assume.

## Citation

Preprint: *From homeostasis to computation: could the GABA switch mature the
core of a signed-XOR motif?* https://zenodo.org/records/21895012.
