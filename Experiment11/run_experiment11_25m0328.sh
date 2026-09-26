#!/usr/bin/env bash
# Experiment 11: computed sulfur-chlorine isomer spectra as the UUV absorber.
# Project roll number: 25m0328
#
# Run from the REPOSITORY ROOT with the venus conda environment active:
#   python3 Experiment11/make_experiment11_inputs_25m0328.py
#   bash Experiment11/run_experiment11_25m0328.sh
#
# Experiment 11 puts the computed photoabsorption cross sections of Trabelsi
# and Francisco (2026) into the radiative transfer of a one-dimensional Venus
# photochemistry transport model. That paper notes that models of this kind
# represent Cl-S chemistry through an empirical opacity rather than explicit
# Cl-S spectra, and does not itself carry its cross sections into such a model.
#
# Protocol matches Experiment 8: spectral shape only, normalised to the same
# band-integrated opacity, nominal vertical placement, A_UUV = 1. The five
# candidate spectra 8A, 8B, 8C, 11A and 11B are therefore directly comparable.

set -euo pipefail

CFG="vulcan_cfg.py"
BACKUP="vulcan_cfg-25m0328-pre-Experiment8.py"
BASELINE="output/Experiment8_baseline_A1-25m0328.vul"

SPECTRA=(
  "atm/UV_absorber_Experiment11A_SSCl2_25m0328.txt"
  "atm/UV_absorber_Experiment11B_ClSSCl_25m0328.txt"
)

# ---------------------------------------------------------------- preflight --
if [ ! -f "$CFG" ]; then
  echo "ERROR: vulcan_cfg.py not found. Run this from the repository root." >&2
  exit 1
fi

missing=0
for f in "${SPECTRA[@]}"; do
  [ -f "$f" ] || { echo "ERROR: missing input $f" >&2; missing=1; }
done
if [ "$missing" -ne 0 ]; then
  echo >&2
  echo "Generate the inputs first:" >&2
  echo "  python3 Experiment11/make_experiment11_inputs_25m0328.py" >&2
  exit 1
fi

if [ ! -f "$BASELINE" ]; then
  echo "ERROR: same-machine baseline $BASELINE not found." >&2
  echo "Experiment 11 is compared against the Experiment 8 baseline. Run that first." >&2
  exit 1
fi

if [ -f "$BACKUP" ]; then
  echo "existing backup kept: $BACKUP"
else
  cp "$CFG" "$BACKUP"
  echo "active config archived as $BACKUP"
fi

restore_cfg () {
  local status=$?
  [ -f "$BACKUP" ] && { cp "$BACKUP" "$CFG"; echo; echo "config restored from $BACKUP"; }
  [ "$status" -ne 0 ] && echo "run stopped with status $status" >&2
  return $status
}
trap restore_cfg EXIT INT TERM

set_cfg () {
  python3 - "$1" "$2" <<'PY'
import re, sys
spectrum, out = sys.argv[1], sys.argv[2]
p = "vulcan_cfg.py"
s = open(p).read()
s = re.sub(r"^N_particle_path = .*$",  "N_particle_path = 'atm/mode1+2.txt'", s, flags=re.M)
s = re.sub(r"^UV_absorber_path = .*$", f"UV_absorber_path = '{spectrum}'",    s, flags=re.M)
s = re.sub(r"^A_UUV = .*$",            "A_UUV = 1.0",                        s, flags=re.M)
s = re.sub(r"^out_name = .*$",         f"out_name =  '{out}'",               s, flags=re.M)
open(p, "w").write(s)
PY
}

run_case () {
  local label="$1" spectrum="$2" out="$3"
  if [ -f "output/$out" ]; then
    echo
    echo "SKIPPING $label: output/$out already exists. Delete it to force a re-run."
    return 0
  fi
  echo
  echo "=============================================================="
  echo " $label"
  echo "=============================================================="
  set_cfg "$spectrum" "$out"
  grep -E "^(N_particle_path|UV_absorber_path|A_UUV|out_name|use_photo|use_condense)" "$CFG"
  echo "--- pre-run gate above, starting integration ---"
  python3 vulcan.py -n
}

run_case "11A: SSCl2 computed spectrum, nominal placement" \
  "atm/UV_absorber_Experiment11A_SSCl2_25m0328.txt" \
  "Experiment11A_SSCl2-25m0328.vul"

run_case "11B: ClSSCl computed spectrum, nominal placement" \
  "atm/UV_absorber_Experiment11B_ClSSCl_25m0328.txt" \
  "Experiment11B_ClSSCl-25m0328.vul"

echo
echo "both runs complete. Next:"
echo "  python3 Experiment11/analyze_experiment11_25m0328.py"
