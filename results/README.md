# Curated results

This directory contains compact outputs from the thesis experiments so the
main conclusions can be inspected without downloading multi-gigabyte snapshot
matrices.

## Tables

`tables/dipole/` contains method errors, operator comparisons, enstrophy peaks,
and online timing summaries. `tables/rayleigh_benard/` contains grid,
method, matched-rank, retained-energy, operator, and prediction comparisons.

The recorded timings are machine- and implementation-dependent. Online
speedups exclude the cost of generating full-order snapshots and building the
reduced model.

## Figures

The selected figures show:

- localized dipole-wall vorticity errors;
- the Rayleigh-Benard circulation reversal;
- equality of exact-OpInf and intrusive trajectories;
- the retained-energy versus online-cost trade-off;
- continuation beyond the ROM training horizon.

The two benchmark images are reference assets used by
`analysis/dipole_validation_plots.jl`.

