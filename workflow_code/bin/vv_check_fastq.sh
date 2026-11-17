#!/usr/bin/env bash
# One-pass gzip CRC + FASTQ record check. Always exit 0.
# stdout: basename<TAB>OK|FAIL<TAB>OK|FAIL<TAB>msg
set -u
src=$1
base=$(basename "$src")

if [[ ! -e "$src" ]]; then
  printf '%s\tFAIL\tFAIL\tmissing\n' "$base"
  exit 0
fi

gzip_err=$(mktemp)
fmt_err=$(mktemp)
set +e
set +o pipefail
gzip -dc "$src" 2>"$gzip_err" | awk '
{
    n++
    r = (n - 1) % 4
    if (!bad && r == 0 && substr($0, 1, 1) != "@") {
        print "header line " n " expected @" > "/dev/stderr"
        bad = 1
    }
    if (!bad && r == 2 && substr($0, 1, 1) != "+") {
        print "sep line " n " expected +" > "/dev/stderr"
        bad = 1
    }
}
END {
    if (n == 0) {
        print "empty FASTQ" > "/dev/stderr"
        bad = 1
    }
    if (n % 4 != 0) {
        print "incomplete record, " n " lines" > "/dev/stderr"
        bad = 1
    }
    if (bad) exit 2
}
' 2>"$fmt_err"
ps=("${PIPESTATUS[@]}")
g_rc=${ps[0]:-1}
a_rc=${ps[1]:-1}
set -e

gzip_st=OK
fmt_st=OK
msg=""

if [[ "$a_rc" -ne 0 ]]; then
  fmt_st=FAIL
  msg=$(tr '\n' ' ' < "$fmt_err")
fi

if [[ "$g_rc" -ne 0 ]]; then
  gzip_st=FAIL
  gmsg=$(tr '\n' ' ' < "$gzip_err")
  msg=$(printf '%s %s' "$msg" "$gmsg")
elif [[ -s "$gzip_err" ]]; then
  gzip_st=FAIL
  gmsg=$(tr '\n' ' ' < "$gzip_err")
  msg=$(printf '%s %s' "$msg" "$gmsg")
fi

msg=$(printf '%s' "$msg" | tr '\t' ' ' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
printf '%s\t%s\t%s\t%s\n' "$base" "$gzip_st" "$fmt_st" "$msg"
rm -f "$gzip_err" "$fmt_err"
exit 0
