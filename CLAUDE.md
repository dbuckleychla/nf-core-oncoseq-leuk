# Project Instructions for AI Agents

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

This file provides instructions and context for AI coding agents working on this project.

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

## Build & Test

_Add your build and test commands here_

```bash
# Example:
# npm install
# npm test
```

## Architecture Overview

_Add a brief overview of your project architecture_

## Conventions & Patterns

_Add your project-specific conventions here_
