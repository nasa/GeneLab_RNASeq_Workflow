#!/usr/bin/env bash
# First N FASTQ records → gzip OUT.
# gzip -dcf: inflate .gz, pass plain FASTQ through.
# 0 = SRC smaller than N. 141 = SIGPIPE after head closed. Both OK if OUT has records.
# Record check: complete quartets, @ header, + separator.
# https://github.com/nextflow-io/nextflow/blob/master/modules/nextflow/src/main/groovy/nextflow/splitter/FastqSplitter.groovy
set -u
nrec=$1
out=$2
src=$3
nlines=$((nrec * 4))

set +e
set -o pipefail
if [[ "$src" == *://* ]]; then
    fetch_uri.sh "$src" - | gzip -dcf | head -n "$nlines" | gzip > "$out"
else
    gzip -dcf "$src" | head -n "$nlines" | gzip > "$out"
fi
rc=$?
set -e

if [[ "$rc" -ne 0 && "$rc" -ne 141 ]]; then
    echo "ERROR: stream failed (exit ${rc}) from ${src}" >&2
    exit 1
fi

# gzip of zero FASTQ lines is still a 20-byte .gz — check inflated content
set -o pipefail
gzip -dcf "$out" | awk -v src="$src" '
{
    n++
    r = (n - 1) % 4
    if (r == 0 && substr($0, 1, 1) != "@") {
        print "ERROR: Invalid FASTQ format (line " n " expected @) from " src > "/dev/stderr"
        bad = 1
        exit 1
    }
    if (r == 2 && substr($0, 1, 1) != "+") {
        print "ERROR: Invalid FASTQ format (line " n " expected +) from " src > "/dev/stderr"
        bad = 1
        exit 1
    }
}
END {
    if (bad) exit 1
    if (n == 0) {
        print "ERROR: empty FASTQ after first records from " src > "/dev/stderr"
        exit 1
    }
    if (n % 4 != 0) {
        print "ERROR: Invalid FASTQ format (incomplete record, " n " lines) from " src > "/dev/stderr"
        exit 1
    }
}
'
