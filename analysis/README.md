# Analysis scripts

These scripts generated the thesis statistics and curated figures. Most expect
full-order and reduced-order run directories under `output/`; those large
binary datasets are not stored in Git.

- `dipole_validation_plots.jl`: grid and benchmark validation figures;
- `dipole_chapter_stats.jl`: principal dipole error tables;
- `dipole_timing_stats.jl`: repeated online-advance timing summary;
- `dipole_light_plots.jl`: compact ROM diagnostic plots;
- `rayleigh_snapshot_plots.jl`: circulation-reversal snapshots;
- `rayleigh_chapter_stats.jl`: switched-ROM comparisons and prediction metrics;
- `compute_rayleigh_paper_comparison.jl` and
  `rayleigh_benchmark_comparison.Rmd`: comparison with published convection
  benchmarks.

The original run identifiers are declared near the start of each script.
Replace them when regenerating experiments with different configuration hashes.
Tables are written to `results/tables/` and figures to `results/figures/`.

