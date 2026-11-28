# Public repository content audit

Audited on 21 November 2026. No access credentials or recognizable real clinical
sample identifiers were found in the inspected public repository content.
Explicitly synthetic fixtures and public reference names remain.

## Scope

- Current tracked files, including the compressed synthetic VCF/BAM fixtures.
- All reachable fork Git references: 2,544 historical file blobs were examined
  for institutional clinical-accession and credential patterns. Commit messages
  were also checked for recognizable clinical accessions.
- Gitleaks 8.30.1 scanned 1,066 reachable Git commits (about 255 MB of patches)
  and the working directory, including archive inspection (about 1.50 GB),
  including initialized submodule source files. The institutional-accession
  guard separately inspected 245 initialized submodule tracked files.
  Archive detection misclassified two plain-text `tab.gz.config` test files;
  those were inspected as text and included in the Git-history scan.
- Public Beads issue text was checked for recognizable clinical accessions.
- Samplesheet examples and generated previews were reviewed: identifiers are
  synthetic, and clinical execution reports/real inputs are not tracked.

Four initial Gitleaks alerts were false positives: three synthetic Secrets
Manager ARNs in Terraform tests, and an image SHA-256 in the synthetic review
preview. `.gitleaks.toml` retains the built-in rules and narrowly excludes those
exact fixture lines; it does not exclude arbitrary credentials, reference
folders or reports. S3 reference locations, reference releases, container/image
hashes, public author attribution and non-secret infrastructure identifiers
are permitted metadata. GitHub Actions secret _references_ contain no values.

## Prevention

`AGENTS.md` and `CLAUDE.md` explicitly prohibit committing or pushing PHI, real
clinical identifiers or actual credentials, including in commit messages,
filenames, issues, reports and fixtures. Local run inputs, output reports and
credential files have additional ignore patterns.

The pre-commit configuration runs a staged clinical-accession/credential guard
and Gitleaks with redacted diagnostics. Install Gitleaks and the repository's
pre-commit runner before contributing. For example, with Homebrew and prek:

```bash
brew install gitleaks
prek install
prek run --all-files
```

Before pushing, inspect the staged changes and messages and run:

```bash
python3 scripts/check_public_content.py --staged
python3 scripts/check_outgoing_secrets.py
```

The outgoing scan checks every new commit patch, not only the final diff; a
secret introduced and then removed still blocks the push. Missing upstream
tracking falls back to full history. Use `--full-history` for an explicit audit,
tag/multi-branch pushes or a changed trust baseline. This avoids repeatedly
scanning the already-audited thousand-commit history during routine fixes.

Do not bypass the checks or use broad allowlists to conceal findings. The guard
reports blob hashes and categories rather than sensitive values or filenames.
It accepts reference metadata and synthetic identifiers; it rejects the known
institutional clinical-accession pattern, including encoded forms and gzip
content. Gitleaks adds broader secret/entropy detection.

Automated scanners cannot establish that arbitrary free text, images, genomic
records or unfamiliar identifiers are free of PHI. Manually verify that every
fixture/preview is synthetic or an approved public reference. Never copy actual
run data into this checkout for publication. No Git-history rewrite was needed
for the inspected findings. A real exposure would require separate history
remediation and credential rotation, where applicable.
