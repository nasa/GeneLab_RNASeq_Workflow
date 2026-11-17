#!/usr/bin/env bash
# src (string) → dest file, or dest '-' for stdout. s3:// uses aws.
set -euo pipefail

insecure=0
if [[ "${1:-}" == "--insecure" ]]; then
    insecure=1
    shift
fi
src=$1
dest=$2

if [[ "$src" == s3://* ]]; then
    aws s3 cp "$src" "$dest" --no-sign-request || aws s3 cp "$src" "$dest"
else
    if [[ "$insecure" -eq 1 ]]; then
        wget --no-check-certificate -q -O "$dest" "$src"
    else
        wget -q -O "$dest" "$src"
    fi
fi

if [[ "$dest" != "-" && ! -s "$dest" ]]; then
    echo "ERROR: empty or missing fetch from ${src}" >&2
    exit 1
fi
