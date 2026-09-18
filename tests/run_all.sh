#!/usr/bin/env bash
# Run all three ground-truth validations against the SHIPPED analysis code and
# write VALIDATION.md. Exits non-zero if any arm fails its thresholds.
# Usage:  bash tests/run_all.sh   (or: make validate)
set -uo pipefail
cd "$(dirname "$0")/.."           # repo root
REPO="$PWD"
echo "== Project 2 validation =="

Rscript tests/validate_de.R;  A1=$?
python3 arm2_singlecell/notebooks/validate_arm2.py >/tmp/arm2.log 2>&1; grep -E '^>>' /tmp/arm2.log | tail -2
Rscript tests/validate_dtu.R; A3=$?

python3 - "$REPO" "$A1" "$A3" <<'PY'
import json, sys, os
repo, a1, a3 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
def load(p):
    try:  return json.load(open(p))
    except Exception: return {}
m1 = load(f"{repo}/tests/.work/arm1/metrics.json")
m3 = load(f"{repo}/tests/.work/arm3/metrics.json")
m2 = load(f"{repo}/results/arm2/validation/metrics.json")
a2_pass = (m2.get("ari",0) >= 0.90) and (m2.get("accuracy",0) >= 0.90)
p1 = (a1 == 0); p3 = (a3 == 0)
rows = [
 ("Arm 1 — bulk DE (DESeq2+edgeR)",
  f"precision {m1.get('precision','?')}, recall {m1.get('recall','?')}, sign {m1.get('sign','?')}",
  "precision>=0.90, recall>=0.80, sign>=0.99", p1),
 ("Arm 2 — scRNA-seq (Leiden+markers)",
  f"ARI {m2.get('ari','?'):.3f}, NMI {m2.get('nmi','?'):.3f}, accuracy {m2.get('accuracy','?'):.3f}" if m2 else "?",
  "ARI>=0.90, accuracy>=0.90", a2_pass),
 ("Arm 3 — long-read DTU (DRIMSeq+stageR)",
  f"precision {m3.get('precision','?')}, recall {m3.get('recall','?')}",
  "precision>=0.80, recall>=0.90", p3),
]
allpass = all(r[3] for r in rows)
with open(f"{repo}/VALIDATION.md","w") as f:
    f.write("# Validation\n\n")
    f.write("Each arm's analysis code is exercised on a simulated dataset with known\n")
    f.write("ground truth by a script under `tests/` (Arm 2's lives with its module).\n")
    f.write("Regenerate this file any time with `make validate`.\n\n")
    f.write("| Arm | Result vs ground truth | Thresholds | Status |\n")
    f.write("|-----|------------------------|------------|--------|\n")
    for name, res, thr, ok in rows:
        f.write(f"| {name} | {res} | {thr} | {'PASS' if ok else 'FAIL'} |\n")
    f.write(f"\n**Overall: {'PASS' if allpass else 'FAIL'}**\n")
print("\n".join(f"  {r[0]:42s} {'PASS' if r[3] else 'FAIL'}" for r in rows))
print("  " + "-"*54)
print(f"  {'OVERALL':42s} {'PASS' if allpass else 'FAIL'}")
sys.exit(0 if allpass else 1)
PY
