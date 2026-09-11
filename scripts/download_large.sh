#!/usr/bin/env bash
# Fetch the scored large instance tier from the GitHub Release.
#
# Usage:
#   scripts/download_large.sh            # the scored large tier
#
# Files arrive gzipped and are decompressed into instances/, so every path in
# the README works unchanged.  Each file's SHA-256 is checked against
# checksums/<tier>.sha256; a mismatch aborts rather than leaving you to debug a
# truncated download as an algorithm bug.

set -euo pipefail

REPO="${RELEASE_REPO:-ythuang0522/shortest-path-competition}"
TAG="${RELEASE_TAG:-v2.0}"

cd "$(dirname "$0")/.."
mkdir -p instances

TIER="${1:-large}"
case "$TIER" in
  large)  tiers=(large) ;;
  *) echo "usage: $0 [large]" >&2; exit 2 ;;
esac

# The scored set. Each instance ships .graph, .queries and .answers -- the
# answer key is distributed with the data because at these query counts nobody
# can afford to regenerate ground truth by running the foundation.
large_names=(road2d_large lattice3d_large local2d_large scalefree_large
             hugeq_large wide64_large)

sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
  else shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

for tier in "${tiers[@]}"; do
  eval "names=(\"\${${tier}_names[@]}\")"
  sums="checksums/${tier}.sha256"
  echo "Downloading ${tier} tier from ${REPO} @ ${TAG} ..."

  for name in "${names[@]}"; do
    for ext in graph queries answers; do
      dest="instances/${name}.${ext}"
      if [[ -f "$dest" ]]; then
        echo "  [skip] ${dest} already present"
        continue
      fi
      url="https://github.com/${REPO}/releases/download/${TAG}/${name}.${ext}.gz"
      echo "  [get ] ${name}.${ext}"
      curl -fL --retry 3 -o "${dest}.gz" "$url"
      gunzip -f "${dest}.gz"

      if [[ -f "$sums" ]]; then
        want=$(awk -v f="${name}.${ext}" '$2 == f {print $1}' "$sums")
        if [[ -n "$want" ]]; then
          got=$(sha256_of "$dest")
          if [[ "$got" != "$want" ]]; then
            echo "  CHECKSUM MISMATCH for ${dest}" >&2
            echo "    expected $want" >&2
            echo "    got      $got" >&2
            rm -f "$dest"
            exit 1
          fi
        fi
      fi
    done
  done
done
echo "Done."
