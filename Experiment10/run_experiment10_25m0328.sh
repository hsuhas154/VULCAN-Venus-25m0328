#!/usr/bin/env bash
# Experiment 10: universality control for the below-cloud response fingerprint.
# Project roll number: 25m0328
#
# Run from the REPOSITORY ROOT with the venus conda environment active:
#   bash Experiment10/run_experiment10_25m0328.sh
#
# Experiments 7B and 9 both place the absorber below the photochemically
# active region and produce sulfur photolysis fingerprints that differ by a
# single scalar (1.0475, spread 0.2% across eight channels), despite having
# different spectra and different vertical distributions.
#
# Two points do not establish a pattern. This experiment adds a third: the S3
# spectrum on the ClS2 vertical profile. If the scaling holds, the below-cloud
# response is a universal fingerprint whose shape is set by the atmosphere and
# whose amplitude alone depends on the absorber. If it does not hold, the
# agreement between 7B and 9 was a coincidence and must not be reported as a
# result.
#
# Both input files already exist in the repository from Experiments 8 and 9.
# No new inputs are generated.

set -euo pipefail

CFG="vulcan_cfg.py"
BACKUP="vulcan_cfg-25m0328-pre-Experiment8.py"

SPECTRUM="atm/UV_absorber_Experiment8C_S3_25m0328.txt"
PROFILE="atm/mode1+2_Experiment9_ClS2profile_25m0328.txt"
OUT="Experiment10_S3-on-ClS2profile-25m0328.vul"
BASELINE="output/Experiment8_baseline_A1-25m0328.vul"

# ---------------------------------------------------------------- preflight --
if [ ! -f "$CFG" ]; then
  echo "ERROR: vulcan_cfg.py not found. Run this from the repository root:" >&2
  echo "  bash Experiment10/run_experiment10_25m0328.sh" >&2
  exit 1
fi

missing=0
for f in "$SPECTRUM" "$PROFILE"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing input file $f" >&2
    missing=1
  fi
done
if [ "$missing" -ne 0 ]; then
  echo >&2
  echo "Both inputs come from Experiments 8 and 9 and should already be in atm/." >&2
  echo "If they are absent, regenerate them with:" >&2
  echo "  python3 Experiment8-9/make_experiment8_9_inputs_25m0328.py" >&2
  exit 1
fi

if [ ! -f "$BASELINE" ]; then
  echo "ERROR: same-machine baseline $BASELINE not found." >&2
  echo "Experiment 10 must be compared against the baseline produced by the" >&2
  echo "Experiment 8 and 9 run. Run that script first." >&2
  exit 1
fi

if [ -f "output/$OUT" ]; then
  echo "output/$OUT already exists. Delete it to force a re-run."
  exit 0
fi

if [ -f "$BACKUP" ]; then
  echo "existing backup kept: $BACKUP"
else
  cp "$CFG" "$BACKUP"
  echo "active config archived as $BACKUP"
fi

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
python3 - "$PROFILE" "$SPECTRUM" "$OUT" <<'PY'
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

echo
echo "=============================================================="
echo " 10: S3 spectrum on the ClS2 vertical distribution"
echo "=============================================================="
grep -E "^(N_particle_path|UV_absorber_path|A_UUV|out_name|use_photo|use_condense)" "$CFG"
echo "--- pre-run gate above, starting integration ---"
python3 vulcan.py -n

echo
echo "run complete. Next:"
echo "  python3 Experiment10/analyze_experiment10_scaling_25m0328.py"
