#!/bin/sh
# Regression check for bin/rcount against the verified values r(1..16), with both
# canonical forms and in --passes mode. Usage: scripts/check_rcount.sh [path/to/rcount]
set -e
RCOUNT=${1:-./bin/rcount}
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
expected="1 2 5 11 33 74 144 232 639 1406 3164 4992 12501 26973 55937 104169"
n=0
for r in $expected; do
  n=$((n + 1))
  for c in nauty library; do
    got=$("$RCOUNT" "$n" --canon "$c" --tmp "$TMP/x" 2>/dev/null | cut -d' ' -f3)
    if [ "$got" != "$r" ]; then echo "FAIL n=$n canon=$c: got $got, expected $r"; exit 1; fi
  done
done
got=$("$RCOUNT" 16 --passes 5 --threads 2 --split 9 --tmp "$TMP/y" 2>/dev/null | cut -d' ' -f3)
if [ "$got" != "104169" ]; then echo "FAIL --passes mode: got $got"; exit 1; fi
echo "rcount: r(1..16) correct with both canonical forms and in --passes mode"
