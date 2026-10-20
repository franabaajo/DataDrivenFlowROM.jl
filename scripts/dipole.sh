#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

JULIA="${JULIA:-}"
EXAMPLES_PROJECT="${EXAMPLES_PROJECT:-$ROOT_DIR}"
DIPOLE_TASKS="${DIPOLE_TASKS:-$ROOT_DIR/experiments/DipoleWallCollision/scripts/dipole_tasks.jl}"
GRID_STUDY="${GRID_STUDY:-$ROOT_DIR/experiments/DipoleWallCollision/configs/grid_study.toml}"
FOM_BASE="${FOM_BASE:-$ROOT_DIR/experiments/DipoleWallCollision/configs/demo_fom.toml}"
ROM_BASE="${ROM_BASE:-$ROOT_DIR/experiments/DipoleWallCollision/configs/demo_rom.toml}"

usage() {
  cat <<EOF
Usage:
  $(basename "$0") run-all [fom_base.toml] [rom_base.toml]
  $(basename "$0") run-fom [fom_base.toml]
  $(basename "$0") sweep-fom [grid_study.toml] [fom_base.toml]
  $(basename "$0") diagnostics-fom [fom_run_id ...]
  $(basename "$0") postprocess-fom [flags]
  $(basename "$0") postprocess-one <fom_run_id> [flags]
  $(basename "$0") plot-10-modes <fom_run_id>
  $(basename "$0") rerun-fom <fom_run_id> <tag>
  $(basename "$0") rerun-fom-cpu <fom_run_id> [tag]
  $(basename "$0") run-rom <fom_run_id> [rom_base.toml]
  $(basename "$0") rerun-rom <rom_run_id>

Environment overrides:
  JULIA, EXAMPLES_PROJECT, DIPOLE_TASKS, GRID_STUDY, FOM_BASE, ROM_BASE

Postprocess flags:
  --force-all
  --force-grid
  --force-vorticity
  --force-vorticity-zoom
  --force-pressure-zoom
  --force-singular-values
  --force-pod-modes
  --force-movie
EOF
}

resolve_julia() {
  if [[ -n "${JULIA:-}" ]]; then
    return
  fi

  if command -v julia >/dev/null 2>&1; then
    JULIA="$(command -v julia)"
    return
  fi

  for candidate in \
    "$HOME/.juliaup/bin/julia" \
    "$HOME/.juliaup/bin/julialauncher" \
    "/usr/local/bin/julia" \
    "/usr/bin/julia" \
    "/opt/julia/bin/julia" \
    "/snap/bin/julia"; do
    if [[ -x "$candidate" ]]; then
      JULIA="$candidate"
      return
    fi
  done

  echo "ERROR: Julia executable not found." >&2
  echo "Set JULIA=/path/to/julia or add julia to your PATH." >&2
  exit 1
}

cmd="${1:-}"
shift || true

case "$cmd" in
  run-all|run-pipeline|pod)
    resolve_julia
    fom_path="${1:-$FOM_BASE}"
    rom_path="${2:-$ROM_BASE}"
    exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" run-all "$fom_path" "$rom_path"
    ;;
  run-fom)
    resolve_julia
    fom_path="${1:-$FOM_BASE}"
    exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" run-fom "$fom_path"
    ;;
  diagnostics-fom)
    resolve_julia
    exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" diagnostics-fom "$@"
    ;;
  sweep-fom|sweep)
    resolve_julia
    grid_path="${1:-$GRID_STUDY}"
    base_path="${2:-$FOM_BASE}"
    exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" sweep-fom "$grid_path" "$base_path"
    ;;
  postprocess-fom|postprocess)
    resolve_julia
    exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" postprocess-fom "$@"
    ;;
  postprocess-one)
    if [[ $# -lt 1 ]]; then
      usage
      exit 1
    fi
    resolve_julia
    fom_id="$1"
    shift
    exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" postprocess-one "$fom_id" "$@"
    ;;
  plot-10-modes)
    if [[ $# -lt 1 ]]; then
      usage
      exit 1
    fi
    resolve_julia
    fom_id="$1"
    exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" plot-10-modes "$fom_id"
    ;;
  rerun-fom)
    if [[ $# -lt 2 ]]; then
      usage
      exit 1
    fi
    resolve_julia
    fom_id="$1"
    tag="$2"
    exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" rerun-fom "$fom_id" "$tag"
    ;;
  rerun-fom-cpu)
    if [[ $# -lt 1 ]]; then
      usage
      exit 1
    fi
    resolve_julia
    fom_id="$1"
    tag="${2:-cpu}"
    DIPOLE_BACKEND=cpu exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" rerun-fom "$fom_id" "$tag"
    ;;
  run-rom)
    if [[ $# -lt 1 ]]; then
      usage
      exit 1
    fi
    resolve_julia
    fom_id="$1"
    rom_path="${2:-$ROM_BASE}"
    exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" run-rom "$fom_id" "$rom_path"
    ;;
  rerun-rom)
    if [[ $# -lt 1 ]]; then
      usage
      exit 1
    fi
    resolve_julia
    rom_id="$1"
    exec "$JULIA" --project="$EXAMPLES_PROJECT" "$DIPOLE_TASKS" rerun-rom "$rom_id"
    ;;
  ""|-h|--help|help)
    usage
    ;;
  *)
    echo "Unknown command: $cmd" >&2
    usage
    exit 1
    ;;
esac
