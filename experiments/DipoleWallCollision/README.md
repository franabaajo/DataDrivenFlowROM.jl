# Dipole-wall collision

This case studies the normal collision of a vortex dipole with a no-slip wall
at `Re = 2500`. It validates the full-order calculation against reference
enstrophy peaks and compares trajectory OpInf, exact OpInf, and matrix-free
intrusive Galerkin ROMs.

The global POD basis captures bulk kinetic-energy decay but has difficulty with
the thin vorticity layers and secondary vortices created at the wall.

## Run the demo

From the repository root:

```bash
DIPOLE_BACKEND=cpu bash scripts/dipole.sh run-all
```

`configs/demo_fom.toml` and `configs/demo_rom.toml` are small smoke-test
configurations. `configs/fom_base.toml` and `configs/rom_base.toml` retain the
thesis-scale settings.

## Contents

- `configs/`: demo, full-scale, grid-sweep, and space-local POD settings;
- `scripts/`: FOM, ROM, sweep, diagnostic, and post-processing commands;
- `studies/`: grid-convergence and enstrophy-validation studies.

Reusable algorithms live in `src/workflows/`; case-specific implementation
lives in `src/dipole/`.
