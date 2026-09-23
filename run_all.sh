#!/usr/bin/env bash
# Regenerate data, run every demo, run every test. One command, whole repo.
set -e
cd "$(dirname "$0")"
echo "=== generating data ==="
python data/make_data.py
for lvl in level1_sequences level2_alignment level3_ml level4_structure; do
  echo; echo "=== $lvl demo ==="
  python "$lvl/solution.py"
done
echo; echo "=== level5_pipeline demo ==="
python level5_pipeline/pipeline.py
echo; echo "=== tests (reference implementations) ==="
LADDER_MOD=solution python -m pytest -q
