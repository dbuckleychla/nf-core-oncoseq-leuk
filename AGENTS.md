# Repository Guidelines

## Public repository: PHI and credentials

Under no circumstances may PHI, real clinical sample/patient identifiers, or
access keys/credentials be committed or pushed to this public repository.
This applies to code, fixtures, reports, logs, filenames, commit messages,
PR text, and the public Beads/Dolt tracker. Do not copy real run outputs or
samplesheets into tracked files. Use explicitly synthetic fixtures and sample
names only; clinical outputs must remain in private runtime storage.

Reference-data paths, reference names, releases, public reference assets,
container digests, and non-secret infrastructure identifiers are allowed.
References to environment variables or GitHub Actions secrets are allowed;
the actual values are never allowed. Keep credentials in the standard AWS
credential chain, environment or Secrets Manager, never in repository files.

Before every commit/push, review staged content and run the public-content and
Gitleaks checks configured in `.pre-commit-config.yaml`. Before a routine
current-branch push, run `python3 scripts/check_outgoing_secrets.py`: it scans
every outgoing commit, including intermediate changes later removed. Use
`--full-history` for audits, tag/multi-branch pushes, or changed trust baselines;
a missing upstream automatically triggers a full-history scan. Do not rescan
already-audited history on every routine push. Never bypass the checks or
add broad exclusions to hide a finding. Inspect generated/binary assets and
commit/issue text manually as well: automated scans cannot identify every form
of PHI. If exposure is found, stop pushing, remove it from pending changes,
and report the need for credential rotation/history cleanup without repeating
the sensitive value.

## Project Structure & Module Organization

`main.nf` is the pipeline entrypoint. Top-level workflows live in
`workflows/`, reusable orchestration is in `subworkflows/`, and tool wrappers
are split between `modules/local/` and `modules/nf-core/`. Configuration is
under `conf/`, including resource presets and execution profiles such as
`conf/profiles/test.config`. Static inputs and templates live in `assets/`,
helper scripts in `bin/`, and sample data in `test_data/`.

Tests mirror the code layout. Use `tests/modules/.../*.nf.test`,
`tests/subworkflows/.../*.nf.test`, and `workflows/tests/*.nf.test` when adding
or updating coverage.

## Build, Test, and Development Commands

Run a minimal local pipeline test with:

```bash
nextflow run . -profile test,docker --outdir test_results
```

Run the remote-assets profile with:

```bash
nextflow run . -profile test_remotes,apptainer --outdir test_results
```

Run nf-test suites with:

```bash
nf-test test --profile debug,test,docker --verbose
```

Lint the pipeline with:

```bash
nf-core pipelines lint .
```

Inspect merged configuration without launching tasks:

```bash
NXF_OFFLINE=true nextflow config -profile test_remotes
```

## Coding Style & Naming Conventions

Use 4-space indentation in `.nf` and `.config` files. Keep comments short and
focused on workflow intent. Follow existing channel naming patterns such as
`ch_input`, `ch_ref`, and `ch_<source>_for_<target>`. Name tests after the file
they cover, for example `tests/subworkflows/local/variant_calling/phasing.nf.test`.

When adding parameters, update `conf/params.config`,
`nextflow_schema.json`, and relevant docs in `README.md` or `docs/`.

## Testing Guidelines

This repository uses `nf-test`. Add or update tests for workflow behavior,
module interfaces, and regression fixes. Prefer the smallest realistic fixture
from `test_data/`. Before opening a PR, run the affected `nf-test` suite and at
least one `nextflow run` profile that exercises the changed path.

## Commit & Pull Request Guidelines

Recent history favors short, imperative commit messages such as `Updated
whatshap`, `Corrected snp calling...`, or issue-linked prefixes like
`NO-54 added option to index single bam file input`. Keep commits scoped to one
change.

Open PRs against `dev`. Include a concise summary, linked issue, test evidence,
and any config or schema updates. If output structure or reports change, note
the user-facing impact clearly.

<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->

## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**

- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->
