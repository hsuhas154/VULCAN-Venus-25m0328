#!/usr/bin/env bash
# Experiments 8 and 9 for the VULCAN-Venus UUV project.
# Project roll number: 25m0328
#
# Run from the REPOSITORY ROOT with the venus conda environment active:
#   bash Experiment8-9/run_experiment8_9_25m0328.sh
#
# The four input files must sit in the repository's own atm/ directory:
#   cp Experiment8-9/atm/*.txt atm/
#
# Each run takes roughly one minute. The active config is archived once, and
# restored automatically on exit, including on error or interrupt.

set -euo pipefail

CFG="vulcan_cfg.py"
BACKUP="vulcan_cfg-25m0328-pre-Experiment8.py"

SPECTRA=(
  "atm/UV_absorber_Experiment8A_ClS2_25m0328.txt"
  "atm/UV_absorber_Experiment8B_OSSO_S2O2_25m0328.txt"
  "atm/UV_absorber_Experiment8C_S3_25m0328.txt"
  "atm/mode1+2_Experiment9_ClS2profile_25m0328.txt"
)

# ---------------------------------------------------------------- preflight --
if [ ! -f "$CFG" ]; then
  echo "ERROR: vulcan_cfg.py not found. Run this from the repository root:" >&2
  echo "  bash Experiment8-9/run_experiment8_9_25m0328.sh" >&2
  exit 1
fi

missing=0
for f in "${SPECTRA[@]}"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing input file $f" >&2
    missing=1
  fi
done
if [ "$missing" -ne 0 ]; then
  echo >&2
  echo "The Experiment 8 and 9 input files belong in the repository's atm/ directory." >&2
  echo "If you copied the bundle in as Experiment8-9/, run:" >&2
  echo "  cp Experiment8-9/atm/*.txt atm/" >&2
  echo "and then start this script again." >&2
  exit 1
fi

# Archive the active config ONCE. Never overwrite an existing backup: on a
# second run that would replace the good original with a partially modified
# config left behind by an interrupted run.
if [ -f "$BACKUP" ]; then
  echo "existing backup kept: $BACKUP"
else
  cp "$CFG" "$BACKUP"
  echo "active config archived as $BACKUP"
fi

# Restore on any exit path: success, error, or interrupt.
restore_cfg () {
  local status=$?
  if [ -f "$BACKUP" ]; then
    cp "$BACKUP" "$CFG"
    echo
    echo "config restored from $BACKUP"
  fi
  if [ "$status" -ne 0 ]; then
    echo "run stopped with status $status" >&2
  fi
  return $status
}
trap restore_cfg EXIT INT TERM

# ------------------------------------------------------------------- runner --
set_cfg () {
  python3 - "$1" "$2" "$3" <<'PY'
import re, sys
profile, spectrum, out = sys.argv[1], sys.argv[2], sys.argv[3]
p = "vulcan_cfg.py"
s = open(p).read()
s = re.sub(r"^N_particle_path = .*$",  f"N_particle_path = '{profile}'",  s, flags=re.M)
s = re.sub(r"^UV_absorber_path = .*$", f"UV_absorber_path = '{spectrum}'", s, flags=re.M)
s = re.sub(r"^A_UUV = .*$",            "A_UUV = 1.0",                     s, flags=re.M)
s = re.sub(r"^out_name = .*$",         f"out_name =  '{out}'",            s, flags=re.M)
open(p, "w").write(s)
PY
}

run_case () {
  local label="$1" profile="$2" spectrum="$3" out="$4"
  if [ -f "output/$out" ]; then
    echo
    echo "=============================================================="
    echo " SKIPPING $label"
    echo " output/$out already exists. Delete it to force a re-run."
    echo "=============================================================="
    return 0
  fi
  echo
  echo "=============================================================="
  echo " $label"
  echo "=============================================================="
  set_cfg "$profile" "$spectrum" "$out"
  grep -E "^(N_particle_path|UV_absorber_path|A_UUV|out_name|use_photo|use_condense)" "$CFG"
  echo "--- pre-run gate above, starting integration ---"
  python3 vulcan.py -n
}

# Baseline on THIS machine. The candidate runs must be compared against a
# reference produced with the same libraries, not against a reference from
# another machine.
run_case "Baseline: nominal spectrum, nominal placement" \
  "atm/mode1+2.txt" "atm/UV_absorber.txt" "Experiment8_baseline_A1-25m0328.vul"

run_case "8A: ClS2 spectrum, nominal placement" \
  "atm/mode1+2.txt" "atm/UV_absorber_Experiment8A_ClS2_25m0328.txt" \
  "Experiment8A_ClS2-25m0328.vul"

run_case "8B: OSSO spectrum, nominal placement" \
  "atm/mode1+2.txt" "atm/UV_absorber_Experiment8B_OSSO_S2O2_25m0328.txt" \
  "Experiment8B_OSSO_S2O2-25m0328.vul"

run_case "8C: S3 spectrum, nominal placement" \
  "atm/mode1+2.txt" "atm/UV_absorber_Experiment8C_S3_25m0328.txt" \
  "Experiment8C_S3-25m0328.vul"

run_case "9: ClS2 spectrum and ClS2 vertical distribution" \
  "atm/mode1+2_Experiment9_ClS2profile_25m0328.txt" \
  "atm/UV_absorber_Experiment8A_ClS2_25m0328.txt" \
  "Experiment9_ClS2-selfconsistent-25m0328.vul"

echo
echo "all five runs complete. Next:"
echo "  python3 Experiment8-9/analyze_experiment8_9_25m0328.py"
echo "  python3 Experiment8-9/candidate_opacity_screen_25m0328.py"
