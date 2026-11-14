"""
group_ln.py
Combine library models from several flies into one group model.

    python group_ln.py --members contrast_blocks/A,contrast_blocks/B --name group/mi1_control \
        [--method reference_sd] [--reference_state 1] [--model_library DIR]

Members must be SwitchedLN models on the same label with the same state values
(or all plain LN). Filters are Baccus & Meister normalised already, so they are
averaged on the first member's lag grid; the pairwise agreement is recorded.
Nonlinearities differ across flies in dF/F scale and offset for reasons that
are not biology (expression, detector), so each member's nonlinearity is first
put on a common footing and then averaged:

  reference_sd   y_m(g) = a_m N(g) + b_m with a_m = the SD and b_m = the mean of
                 member m's binned response in the REFERENCE state (the high
                 contrast by default). One (a_m, b_m) per fly, shared by all its
                 states, so each fly's gain ratio between states survives. The
                 group nonlinearity is mean(a) N + mean(b), in average dF/F;
                 per-member a and b are stored so a prediction can be put back
                 in any one fly's units.
  plain          average the nonlinearities as they are.

The group is formed where every member has data: the x grid spans the
intersection of the members' measured ranges in each state.
"""
import argparse
import datetime

import numpy as np

from visanalysis.model import spec as mspec
from visanalysis.model.models import model_from_spec, SwitchedLN, LN


def _states(m):
    return dict(m.states) if isinstance(m, SwitchedLN) else {"all": m}


def group_models(specs, name, method="reference_sd", reference_state=1, n_grid=24):
    models = [model_from_spec(s) for s in specs]
    S = [_states(m) for m in models]
    keys = list(S[0])
    for s in S[1:]:
        if list(s) != keys:
            raise ValueError("members have different states: {} vs {}".format(keys, list(s)))
    ref = reference_state if reference_state in keys else keys[-1]
    # per-member scale and offset from the reference state's total nonlinearity
    ab = []
    for s in S:
        x, y, _ = s[ref].nonlinearity.total()
        ab.append((float(y.std()) if method == "reference_sd" else 1.0, float(y.mean()) if method == "reference_sd" else 0.0))
    a_bar, b_bar = float(np.mean([a for a, _ in ab])), float(np.mean([b for _, b in ab]))
    states, report = [], {}
    for v in keys:
        # filter: average densities on the first member's grid
        f0 = S[0][v].filter
        edges = np.asarray(f0.lag_edges) if f0.kind == "tabulated" else np.arange(int(f0.length / 0.01) + 1) * 0.01
        lag_dt = float(np.min(np.diff(edges)))
        dens = []
        for s in S:
            f = s[v].filter
            _, d = f.tabulate(lag_dt)
            d = np.r_[d, np.zeros(max(edges.size - 1 - d.size, 0))][:edges.size - 1]
            dens.append(d)
        dens = np.array(dens)
        corr = [float(np.corrcoef(dens[i], dens[j])[0, 1]) for i in range(len(S)) for j in range(i + 1, len(S))]
        h = dens.mean(0)
        # nonlinearity: normalise each member, average on a common grid inside every member's range
        lo = max(s[v].nonlinearity.range[0] for s in S)
        hi = min(s[v].nonlinearity.range[1] for s in S)
        xg = np.linspace(lo, hi, n_grid)
        N = []
        for s, (a, b) in zip(S, ab):
            x, y, _ = s[v].nonlinearity.total()
            N.append((np.interp(xg, x, y) - b) / a)
        N = np.array(N)
        yg = a_bar * N.mean(0) + b_bar
        se = a_bar * N.std(0) / np.sqrt(len(S))
        states.append(dict(value=v if v != "all" else None, name=S[0][v].name,
                           filter=dict(kind="tabulated", lag_edges=edges, density=h, output_units=f0.output_units,
                                       normalization=f0.normalization, name=S[0][v].name),
                           nonlinearity=dict(kind="tabulated", x=xg, y=yg, y_se=se, offset=0.0,
                                             input_units=S[0][v].nonlinearity.input_units, extrapolation="linear",
                                             name=S[0][v].name)))
        report[str(v)] = dict(filter_pairwise_r=corr, nl_range=[float(lo), float(hi)],
                              filter_variance_ratio=float((h ** 2).sum() / np.mean((dens ** 2).sum(1))))
    m0 = models[0]
    desc = "Group of {} models ({}): mean normalised filter; nonlinearities combined by '{}' on state {}".format(
        len(models), ", ".join(s["name"] for s in specs), method, ref)
    out = mspec.new_spec(name=name.rsplit("/", 1)[-1], kind=m0.kind, source="group", description=desc,
                         channel=m0.output.get("channel"), signal=m0.output.get("signal", "imaged"),
                         units=m0.output.get("units", "dF/F"), indicator=m0.output.get("indicator"),
                         method=method, reference_state=ref,
                         members=[dict(name=s.get("name"), path=s.get("_path"), hash=mspec.model_hash(s),
                                       a=a, b=b, source=s.get("provenance", {}).get("source"),
                                       fly=s.get("provenance", {}).get("fly"))
                                  for s, (a, b) in zip(specs, ab)],
                         a_mean=a_bar, b_mean=b_bar, per_state=report)
    if m0.kind == "SwitchedLN":
        out["label"] = m0.label
        out["states"] = states
        out["default"] = None
    else:
        out["filter"] = states[0]["filter"]
        out["nonlinearity"] = states[0]["nonlinearity"]
    return out, report


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--members", required=True, help="comma-separated library names")
    ap.add_argument("--name", required=True, help="library name for the group model, e.g. group/mi1_control")
    ap.add_argument("--method", default="reference_sd", choices=("reference_sd", "plain"))
    ap.add_argument("--reference_state", type=int, default=1)
    ap.add_argument("--model_library", default=None)
    ap.add_argument("--overwrite", default="True")
    args = ap.parse_args()
    specs = [mspec.load_model(n.strip(), library=args.model_library) for n in args.members.split(",") if n.strip()]
    out, report = group_models(specs, args.name, args.method, args.reference_state)
    p = mspec.save_model(out, args.name, library=args.model_library, overwrite=args.overwrite == "True")
    print("group model {} from {} members -> {}".format(args.name, len(specs), p))
    for m in out["provenance"]["members"]:
        print("  {}: a {:.4f}, b {:+.4f}".format(m["name"], m["a"], m["b"]))
    for v, r in report.items():
        print("  state {}: filter pairwise r {}, NL range [{:.3f}, {:.3f}], filter variance ratio {:.3f}".format(
            v, ", ".join("{:.3f}".format(c) for c in r["filter_pairwise_r"]), *r["nl_range"], r["filter_variance_ratio"]))


if __name__ == "__main__":
    main()
