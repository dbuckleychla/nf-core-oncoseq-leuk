#!/usr/bin/env bash
# Replay one Batch basecalling shard locally. No result uploads or AWS writes.
set -euo pipefail
umask 077

PLAN_ONLY=false
if [[ ${1:-} == --plan-only ]]; then PLAN_ONLY=true; shift; fi
if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 [--plan-only] BATCH_JOB_ID [LOCAL_WORK_DIRECTORY]" >&2
    exit 2
fi
JOB_ID=$1
WORK_DIR=${2:-"$HOME/oncoseq-shard-replay/$JOB_ID"}
AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
export AWS_DEFAULT_REGION
for tool in aws python3; do command -v "$tool" >/dev/null; done
mkdir -p "$WORK_DIR"
WORK_DIR=$(cd "$WORK_DIR" && pwd -P)
if git -C "$WORK_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo 'Use a work directory outside any Git checkout; shard files may contain PHI.' >&2
    exit 1
fi
if [[ -e "$WORK_DIR/.command.sh" ]]; then
    echo 'Use a fresh work directory to avoid overwriting replay inputs or outputs.' >&2
    exit 1
fi

aws batch describe-jobs --jobs "$JOB_ID" --output json > "$WORK_DIR/job.json"
python3 - "$WORK_DIR" <<'PY'
import json, re, sys
from pathlib import Path
root = Path(sys.argv[1])
jobs = json.loads((root/'job.json').read_text()).get('jobs', [])
if len(jobs) != 1:
    raise SystemExit('Batch job not found; check region and AWS profile')
container = jobs[0]['container']
command = '\n'.join(container.get('command', []))
match = re.search(r's3://[^\s"\x27;]+/\.command\.run', command)
if not match:
    raise SystemExit('Cannot locate Nextflow work wrapper in this job')
image = container.get('image', '')
if not re.fullmatch(r'.+@sha256:[a-f0-9]{64}', image):
    raise SystemExit('Replay requires a digest-pinned image')
resources = {r['type']: r['value'] for r in container.get('resourceRequirements', [])}
cpus = resources.get('VCPU', container.get('vcpus'))
memory = resources.get('MEMORY', container.get('memory'))
if not cpus or not memory or int(resources.get('GPU', 0)) != 1:
    raise SystemExit('Expected a one-GPU Batch job with explicit CPU/memory requests')
(root/'work-uri.txt').write_text(match[0].rsplit('/', 1)[0]+'\n')
(root/'image.txt').write_text(image+'\n')
(root/'limits.txt').write_text(str(cpus)+'\n'+str(memory)+'\n')
PY
WORK_URI=$(cat "$WORK_DIR/work-uri.txt")
aws s3 cp "$WORK_URI/.command.run" "$WORK_DIR/.command.run" --only-show-errors
aws s3 cp "$WORK_URI/.command.sh" "$WORK_DIR/.command.sh" --only-show-errors

# Parse staging descriptions as data. Never execute .command.run: it uploads outputs.
python3 - "$WORK_DIR" <<'PY'
import json, shlex, sys
from pathlib import Path
root = Path(sys.argv[1])
script = (root/'.command.sh').read_text()
if 'dorado basecaller' not in script:
    raise SystemExit('This is not a Dorado basecalling shard')
rows = []
seen = set()
for line in (root/'.command.run').read_text().splitlines():
    line = line.strip()
    if not line.startswith('downloads+=('):
        continue
    quoted = shlex.split(line[len('downloads+=('):-1])
    if len(quoted) != 1:
        raise SystemExit('Unsupported staging syntax')
    parts = shlex.split(quoted[0])
    if len(parts) != 3 or parts[0] != 'nxf_s3_download':
        raise SystemExit('Unsupported staging command')
    source, dest = parts[1:]
    if dest in ('.command.sh', '.command.run'):
        continue
    path = Path(dest)
    if not source.startswith('s3://') or path.is_absolute() or '..' in path.parts:
        raise SystemExit('Unsafe input destination or unsupported source')
    if dest in seen:
        raise SystemExit('Duplicate staging destination')
    seen.add(dest)
    directory = path.suffix not in ('.pod5', '.bam') and path.name not in ('NOFILE', 'NOMOD')
    rows.append({'source': source, 'destination': dest, 'directory': directory})
if not rows or not any(r['destination'].endswith('.pod5') for r in rows):
    raise SystemExit('No POD5 shard staging descriptions found')
(root/'staging.json').write_text(json.dumps(rows, indent=2)+'\n')
with (root/'staging.tsv').open('w') as stream:
    for row in rows:
        if any(c in row['source']+row['destination'] for c in '\t\n'):
            raise SystemExit('Unsupported whitespace in staging paths')
        stream.write(row['source']+'\t'+row['destination']+'\t'+str(int(row['directory']))+'\n')
print('Replay plan:', sum(r['destination'].endswith('.pod5') for r in rows),
        'POD5 files;', sum(r['directory'] for r in rows), 'model directories')
PY
IMAGE=$(cat "$WORK_DIR/image.txt")
CPU_LIMIT=$(sed -n '1p' "$WORK_DIR/limits.txt")
MEMORY_MIB=$(sed -n '2p' "$WORK_DIR/limits.txt")
printf 'Image: %s\nOriginal task limits: %s CPUs / %s MiB / one GPU\n' "$IMAGE" "$CPU_LIMIT" "$MEMORY_MIB"
if [[ $PLAN_ONLY == true ]]; then
    echo 'Plan only: downloaded job metadata/wrapper/command; no POD5/model downloads or execution.'
    exit 0
fi
for tool in docker nvidia-smi; do command -v "$tool" >/dev/null; done

while IFS=$'\t' read -r source destination directory; do
    if [[ $directory == 1 ]]; then
        mkdir -p "$WORK_DIR/$destination"
        aws s3 cp "${source%/}/" "$WORK_DIR/$destination/" --recursive --only-show-errors
    else
        mkdir -p "$WORK_DIR/$(dirname "$destination")"
        aws s3 cp "$source" "$WORK_DIR/$destination" --only-show-errors
    fi
done < "$WORK_DIR/staging.tsv"

docker pull --platform linux/amd64 "$IMAGE"
CONTAINER_NAME="oncoseq-replay-${JOB_ID}-$$"
MONITOR_PIDS=()
cleanup() {
    for pid in "${MONITOR_PIDS[@]}"; do kill "$pid" 2>/dev/null || true; done
    docker stop "$CONTAINER_NAME" >/dev/null 2>&1 || true
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
nvidia-smi --query-gpu=timestamp,index,name,utilization.gpu,utilization.memory,memory.used,memory.total,power.draw,power.limit,clocks.sm,clocks.mem,temperature.gpu \
    --format=csv -l 1 > "$WORK_DIR/gpu.csv" 2> "$WORK_DIR/gpu-monitor.err" &
MONITOR_PIDS+=("$!")
if command -v vmstat >/dev/null; then
    vmstat -t 1 > "$WORK_DIR/vmstat.log" & MONITOR_PIDS+=("$!")
fi
if command -v iostat >/dev/null; then
    iostat -x -t 1 > "$WORK_DIR/iostat.log" & MONITOR_PIDS+=("$!")
fi
printf 'Start: %s\n' "$(date -u +%FT%TZ)" > "$WORK_DIR/replay-timing.txt"
# The container has no network or AWS credentials. Outputs stay in this directory.
docker run --rm --name "$CONTAINER_NAME" --platform linux/amd64 \
    --gpus "device=${REPLAY_GPU:-0}" --network none \
    --cpus "$CPU_LIMIT" --memory "${MEMORY_MIB}m" \
    -v "$WORK_DIR:/work" -w /work --entrypoint bash "$IMAGE" \
    -euo pipefail .command.sh 2>&1 | tee "$WORK_DIR/basecall.log"
printf 'End: %s\n' "$(date -u +%FT%TZ)" >> "$WORK_DIR/replay-timing.txt"
echo "Finished. BAM, versions and monitoring files are local in $WORK_DIR. Nothing was uploaded."
