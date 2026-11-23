#!/usr/bin/env bash
# Inspect only the first 1000 records; never log read names or sequences.
set -euo pipefail
if [[ $# != 3 ]]; then
    echo 'Usage: check_ubam.sh INPUT_BAM OUTPUT_JSON REQUIRED_TAGS_CSV' >&2
    exit 2
fi
ubam_input=$1
ubam_report=$2
ubam_tags=$3
samtools quickcheck -u "$ubam_input"
set +e
samtools view "$ubam_input" | awk -v tags="$ubam_tags" '
BEGIN { FS="\t"; n=split(tags,required,","); checked=0; bad=0 }
{
    checked++
    flag=$2+0
    if (int(flag/4)%2 != 1 || int(flag/256)%2 || int(flag/2048)%2) {
        print "uBAM preflight: sampled record " checked " is not an unmapped primary read" > "/dev/stderr"
        bad=1; exit 1
    }
    if ($10=="*" || $11=="*") {
        print "uBAM preflight: sampled record " checked " lacks sequence or base qualities" > "/dev/stderr"
        bad=1; exit 1
    }
    for (j=1;j<=n;j++) {
        found=0
        for (i=12;i<=NF;i++) if (substr($i,1,3)==required[j] ":") found=1
        if (!found) {
            print "uBAM preflight: sampled record " checked " lacks required tag " required[j] > "/dev/stderr"
            bad=1; exit 1
        }
    }
    if (checked==1000) exit
}
END {
    if (bad) exit 1
    if (!checked) { print "uBAM preflight: input has no reads" > "/dev/stderr"; exit 1 }
    printf "{\"checked_reads\":%d,\"check_limit\":1000,\"required_tags\":[",checked
    for (j=1;j<=n;j++) printf "%s\"%s\"",(j>1 ? "," : ""),required[j]
    print "],\"sampled_only\":true,\"remaining_reads_assumed_tagged\":true}"
}' > "$ubam_report"
ubam_pipe_status=("${PIPESTATUS[@]}")
set -e
# The producer may receive SIGPIPE when the sampler stops at record 1000.
if [[ ${ubam_pipe_status[1]} != 0 || ( ${ubam_pipe_status[0]} != 0 && ${ubam_pipe_status[0]} != 141 ) ]]; then
    rm -f "$ubam_report"
    echo 'uBAM preflight failed; alignment was not started' >&2
    exit 1
fi
