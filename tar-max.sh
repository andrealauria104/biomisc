#!/usr/bin/env bash
# tar-max.sh
#
# Usage: ./tar-max.sh <xz|7z> <folder> <archive-base-name> [threads]
#
# Examples:
#   ./tar-max.sh xz mydir backup_xz 4
#   ./tar-max.sh 7z mydir backup_7z 8
#
# Produces either backup_xz.tar.xz or backup_7z.tar.7z

set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 <xz|7z> <folder> <archive-base-name> [threads]" >&2
  exit 1
fi

METHOD="$1"
FOLDER="$2"
BASENAME="$3"
THREADS="${4:-1}"

case "$METHOD" in
  xz)
    OUT="${BASENAME}.tar.xz"
    tar -I "xz -9e -T${THREADS}" -cf "$OUT" "$FOLDER"
    echo "Done: $OUT"
    echo "Extract with: tar -xJf \"$OUT\""
    ;;
  7z)
    OUT="${BASENAME}.tar.7z"
    tar -cf - "$FOLDER" | 7z a -t7z -mx=9 -mmt="$THREADS" -si "$OUT" >/dev/null
    echo "Done: $OUT"
    echo "Extract with: 7z x \"$OUT\" && tar -xf \"${BASENAME}.tar\""
    ;;
  *)
    echo "Error: method must be 'xz' or '7z'" >&2
    exit 1
    ;;
esac
