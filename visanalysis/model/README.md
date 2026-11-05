# visanalysis.model — forward LN models

Run LN models forward on stimuli to see what an experiment should show, and to
test whether an LN description holds. The analysis pipelines in
`visanalysis/util` *estimate* models from data; this package *runs* them. The
two meet only at the **model library**.

## The pipeline, in one line

```
stimulus (contrast + labels) -> model -> [indicator] -> sampling + noise -> readouts
```

| stage | module | interface |
|---|---|---|
| stimulus | `stimulus.py` | `Stimulus(t, c, labels)` from a protocol's parameters, a recording, or a synthetic signal |
| model | `models.py`, `filters.py`, `nonlinearities.py` | `model.respond(stim) -> Response` |
| measurement | `measurement.py` | `Indicator.apply` (neural models only), `sample(...)` at imaging times with noise |
| readouts | `readouts.py`, `plots.py` | switch-triggered average, per-state stats, `recover_ln` (fit the model back), summary figure |
| library | `spec.py`, `export.py` | YAML + npz per model; pipelines export, predictor loads |

Models key on **labels** (`contrast_state`, `since_switch`, ...) that the
stimulus carries, so nothing in the model code knows what experiment it is in.

## The library

One YAML (+ npz for arrays) per model, in a directory resolved in this order:

1. an explicit path (`--model_library`, or `library=` in code)
2. `$VISANALYSIS_MODEL_LIBRARY`
3. `model_library:` in `~/.visanalysis/config.yaml`

There is no default. Measured, group and hypothetical models share one format;
provenance says which is which. Hand-written examples are in `examples/`.

## Conventions worth knowing

* Filters are stored as **densities** (response per unit contrast per second), on
  the lag intervals they were fitted on. Grid-independent; a pipeline tap filter
  is `density * lag_dt`.
* Baccus & Meister (2002): a filter with `output_units: contrast` was scaled so
  `var(F*s) = var(s)` on its fitting stimulus; the nonlinearity carries the gain.
  The scale is applied once at export. `bm_normalize(model, stim)` does the same
  for a hypothetical model — never for an exported one.
* A filter's `output_units` must equal its nonlinearity's `input_units`.
* `output.signal` is `imaged` (fitted from fluorescence; no indicator to apply) or
  `neural` (an indicator kernel may be applied). Applying an indicator to an
  imaged model is refused.
* A tabulated nonlinearity knows the generator range it was measured on; every
  prediction reports how much of it ran outside that range.
* `SwitchedLN` applies the current state's LN to the *whole* stimulus history,
  so every parameter changes the instant the label does. That is the simplest
  hypothesis about adaptation, and the switch-triggered average is where it fails
  if it fails.

## Typical use

```
# export fitted models while running the contrast-switch pipeline
python analyze_data_trf_contrast_switch.py ... --export_models True

# predict on a new design
python predict_ln.py --models contrast_switch/JS256_x_JS261_fly_011_series_002_calcium \
    --stimulus ContrastSwitch --stim_params "low_contrast=0.2,high_contrast=0.92,modulator_dwell=5,probe_tau=0.2" \
    --duration 300 --noise_sd 0.05 --output_directory out

# predict on what a recording actually showed
python predict_ln.py --models <name> --stimulus recording --experiment_file_directory <fly folder> --series 2 \
    --output_directory out

# a hypothetical model, put in measured units on this stimulus
python predict_ln.py --models visanalysis/model/examples/tm9_like_off_switched.yaml --bm_normalize True ...
```

Tests: `python -m pytest tests/test_model.py -q` (needs `clandinin_labpack` beside this checkout).
