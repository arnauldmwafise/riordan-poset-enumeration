#!/bin/sh
# Regression check for the independent counter against the verified values r(1..14),
# with several threads, in disk mode and in pass mode.
# Usage: scripts/check_independent.sh [path/to/independent_rcount]
set -e
IRC=${1:-./bin/independent_rcount}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
expected="1 2 5 11 33 74 144 232 639 1406 3164 4992 12501 26973"
n=0
for r in $expected; do
  n=$((n + 1))
  got=$("$IRC" "$n" "$TMP/d" --threads 3 --split 6 2>/dev/null | cut -d' ' -f3)
  if [ "$got" != "$r" ]; then echo "FAIL n=$n (disk mode): got $got, expected $r"; exit 1; fi
done
got=$("$IRC" 14 --passes 3 --threads 4 --split 8 2>/dev/null | cut -d' ' -f3)
if [ "$got" != "26973" ]; then echo "FAIL pass mode: got $got"; exit 1; fi
echo "independent_rcount: r(1..14) correct with threads, in disk and pass mode"
