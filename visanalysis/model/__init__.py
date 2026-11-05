"""Forward models of visual responses: stimulus -> model -> measurement -> readout.

This package RUNS models forward. The pipelines in visanalysis.util ESTIMATE
models from data. The two meet only at the model library (spec.py): a pipeline
exports what it fitted, a group script combines exported models, a hypothetical
model is written by hand, and the predictor reads any of them without knowing
which it is.

Stages and their interfaces
---------------------------
stimulus.Stimulus      contrast over time on a fine grid, plus named LABELS
                       (contrast_state, dwell, since_switch, ...) that models
                       may key on. Sources: a protocol's parameters (before any
                       fly is run), a recording (the frames the display showed),
                       or synthetic signals.
models.Model           .respond(stimulus) -> Response. LN is one filter and one
                       nonlinearity; SwitchedLN picks an LN by a label's value.
                       Nothing in here knows what "contrast switch" is.
measurement            an optional indicator kernel (only for a model whose
                       output is NEURAL rather than already AS IMAGED), sampling
                       at imaging times, additive noise.
readouts               switch-triggered averages, per-label statistics, and
                       recovery of the filter and nonlinearity from the predicted
                       trace, which is the check that the machinery is right.

Conventions
-----------
* Contrast is Weber contrast (I - I0) / I0 and is 0 in the grey.
* A filter is stored as a DENSITY: response per unit contrast per second of
  stimulus, on lag intervals. That makes it independent of the grid it was
  fitted or applied on. A pipeline filter "per tap of lag_dt" is density * lag_dt.
* An LN whose filter outputs 'contrast' units is in the Baccus & Meister (2002)
  convention: the filter was scaled so var(F * s) = var(s) on the stimulus it
  was fitted to, so the nonlinearity carries all the gain. The scale is applied
  ONCE, at export, and stored; prediction never re-normalises, so a model gives
  the same response to the same stimulus whatever else it is shown.
* A filter and its nonlinearity are a matched pair: the nonlinearity's input
  units must equal the filter's output units.
"""
