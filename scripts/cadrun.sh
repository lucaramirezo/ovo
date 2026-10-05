#!/usr/bin/env bash
# Run a geometry build / check with at most 2 such jobs alive machine-wide, so parallel
# agents do not stack OCC processes. Usage: scripts/cadrun.sh .venv/bin/python -m ovo.cad.blade
cd "$(dirname "$0")/.." || exit 1
command -v flock > /dev/null || { echo 'cadrun.sh needs flock (util-linux; macOS: brew install flock)' >&2; exit 127; }
mkdir -p out/work
while true; do
  for i in 1 2; do
    exec {fd}>"out/work/cad$i.lock"
    if flock -n "$fd"; then
      "$@"; rc=$?
      flock -u "$fd"; exec {fd}>&-
      exit $rc
    fi
    exec {fd}>&-
  done
  sleep 3
done
