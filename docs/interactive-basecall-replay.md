# Interactive, local-only basecalling shard replay

Use `scripts/replay_basecall_shard.sh` on a Linux NVIDIA instance with Docker,
NVIDIA Container Toolkit, AWS CLI and Python 3. The instance role or usual AWS
profile needs read access to the Batch job, its work objects, POD5 inputs and
model assets. Allocate local scratch for that shard and the container image. Use an idle
GPU and the same instance class for the most comparable performance test.

```bash
bash scripts/replay_basecall_shard.sh BATCH_JOB_ID /scratch/dorado-replay
```

Use a directory outside Git and a fresh directory for each replay. These are
real clinical inputs/results; do not commit them, their names, or their logs.
On an EC2 host with the appropriate instance role, no named AWS profile is
necessary. To use an existing configured read-only profile, set `AWS_PROFILE`
when invoking the helper. `AWS_DEFAULT_REGION` defaults to `us-west-2`.

The helper downloads the job metadata and Nextflow command/staging descriptions,
parses the staged paths as data, and downloads that exact shard's POD5 files,
offline model directories and optional resume input. It uses the job's immutable
container digest and original CPU, memory and one-GPU reservation. `REPLAY_GPU`
selects the host GPU index, defaulting to zero.

It runs only `.command.sh`. **It never executes `.command.run`**, which contains
Nextflow upload/finalisation routines. The replay container has no network and
receives no AWS credentials. All outputs remain local. The helper makes no
Batch changes, sends no SSM commands and uploads no results.

Planning without downloading POD5/model assets or starting Docker is available:

```bash
bash scripts/replay_basecall_shard.sh --plan-only BATCH_JOB_ID /scratch/dorado-plan
```

Planning still downloads small job/wrapper/command metadata to the specified
private directory. Choose a new directory for the actual replay.

During execution, the helper records `gpu.csv`, `vmstat.log` and `iostat.log`
(the last two when the tools are installed), plus `basecall.log` and UTC timing.
GPU monitoring captures actual GPU duty cycle, memory activity/VRAM, power,
clocks and temperature. For a second terminal:

```bash
tail -f /scratch/dorado-replay/gpu.csv
```

To inspect per-thread host CPU usage while Dorado is running:

```bash
replay_pid=$(pgrep -n -f '[d]orado basecaller')
top -H -p "$replay_pid"
```

Compare GPU activity with power, CPU threads, I/O wait and disk activity. High
GPU duty cycle with modest power can be memory-bound or otherwise less
compute-intensive work. Low duty cycle with one busy host thread suggests a
serial feeding/decoding stage. This initial replay preserves batching and model
options; change one parameter at a time only in a separate controlled benchmark.
