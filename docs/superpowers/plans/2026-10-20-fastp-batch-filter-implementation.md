# GEO FASTQ Batch fastp Filtering Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build, test, and deploy one safe Bash script that filters all 66 GEO SRRs with the approved per-GSE fastp rules while preserving every original FASTQ.

**Architecture:** The script reads the existing six-column `srr_manifest.tsv`, validates each dataset and its expected SE/PE layout, and runs one 60-thread fastp process at a time. Each SRR is written to a sibling staging directory and atomically renamed to `fastp/` only after non-empty FASTQ and report checks pass.

**Tech Stack:** Bash 5, fastp 1.0.1, Git for Windows Bash for automated tests, PowerShell/scp/ssh for deployment.

## Global Constraints

- Remote root: `/home/guanguiwen/10T4/murinematchtx_GEO`.
- fastp binary: `/home/guanguiwen/.conda/envs/fastp/bin/fastp`.
- Use `--thread 60 --qualified_quality_phred 20 --unqualified_percent_limit 40 --n_base_limit 5 --length_required 30 --disable_trim_poly_g`.
- Never trim the 5′ end; never deduplicate; never enable low-complexity filtering, polyX trimming, correction, or merging.
- Enable `--detect_adapter_for_pe` only for GSE285835 and GSE118802; disable adapter trimming for the other seven GSEs.
- Apply `--trim_tail1 1 --trim_tail2 1` only to GSE70046.
- Never overwrite, move, rename, or delete input FASTQ files; do not start the 189 GiB production run.

---

### Task 1: Add executable behavior tests

**Files:**
- Create: `tests/test_run_fastp_batch.sh`
- Test target: `run_fastp_batch.sh`

**Interfaces:**
- Consumes: `BASE_DIR`, `FASTP_BIN`, and an optional manifest path passed as argument 1.
- Produces: a zero exit status only when layout, arguments, output isolation, skip behavior, and partial-output protection are correct.

- [ ] **Step 1: Write the failing test**

Create a Bash test that builds four fixture rows using real SRR/GSE identifiers:

```bash
1.COAD  GSE263502  SRR28591489  # SE, no adapter
1.COAD  GSE245724  SRR26425828  # PE, no adapter
1.COAD  GSE285835  SRR31881982  # PE, adapter detection
2.LUAD  GSE70046   SRR2071501   # PE, trim one tail base
```

The test provides a fake fastp executable that records arguments and creates the requested output, JSON, and HTML files. Assertions must require:

```bash
grep -F -- "--thread 60" "$CALL_LOG"
grep -F -- "--qualified_quality_phred 20" "$CALL_LOG"
grep -F -- "--unqualified_percent_limit 40" "$CALL_LOG"
grep -F -- "--length_required 30" "$CALL_LOG"
grep -F -- "GSE285835" "$CALL_LOG"
grep -F -- "--detect_adapter_for_pe" "$CALL_LOG"
grep -F -- "--trim_tail1 1 --trim_tail2 1" "$CALL_LOG"
! grep -E -- "--trim_front|--cut_front|--dedup|--low_complexity_filter" "$CALL_LOG"
```

The test must also compare SHA-256 hashes of all original fixture FASTQs before and after execution, verify `.complete`, verify a completed SRR is skipped, and verify a pre-existing incomplete `fastp/` directory causes non-zero exit.

- [ ] **Step 2: Run the test and verify RED**

Run:

```powershell
& 'E:\Git\Git\bin\bash.exe' tests/test_run_fastp_batch.sh
```

Expected: non-zero exit with `run_fastp_batch.sh is missing`.

### Task 2: Implement the batch script

**Files:**
- Create: `run_fastp_batch.sh`
- Test: `tests/test_run_fastp_batch.sh`

**Interfaces:**
- Consumes: canonical six-column TSV manifest; raw FASTQs named `<SRR>.fastq.gz` or `<SRR>_1.fastq.gz` plus `<SRR>_2.fastq.gz`.
- Produces: `<SRR>/fastp/` containing filtered FASTQ(s), `fastp.json`, `fastp.html`, `fastp.log`, and `.complete`.

- [ ] **Step 1: Implement strict setup and manifest validation**

Use:

```bash
#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${BASE_DIR:-$SCRIPT_DIR}"
FASTP_BIN="${FASTP_BIN:-/home/guanguiwen/.conda/envs/fastp/bin/fastp}"
MANIFEST="${1:-${SCRIPT_DIR}/srr_manifest.tsv}"
THREADS=60
```

Require an executable fastp, readable manifest, valid cancer/GSE/SRR fields, unique SRRs, and exactly 66 rows when the default production manifest is used.

- [ ] **Step 2: Implement dataset policy mapping**

Map each `cancer/gse` to one expected layout and policy:

```bash
1.COAD/GSE263502|2.LUAD/GSE163087) layout=SE ;;
1.COAD/GSE285835|3.BRCA/GSE118802) layout=PE; adapter=detect ;;
2.LUAD/GSE70046) layout=PE; tail_trim=1 ;;
1.COAD/GSE245724|2.LUAD/GSE109629|3.BRCA/GSE159203|3.BRCA/GSE262136) layout=PE ;;
*) exit with an unsupported-dataset error ;;
```

Default `adapter=disable` and `tail_trim=0` for every row before applying the case.

- [ ] **Step 3: Implement isolated fastp execution**

Build a Bash argument array containing the approved shared parameters and `--dont_overwrite`. Add `--detect_adapter_for_pe` or `--disable_adapter_trimming` from policy; add GSE70046 tail options only for that dataset. Create a staging directory with `mktemp -d`, redirect fastp output to `fastp.log`, require all outputs and reports to be non-empty, create `.complete`, then rename the staging directory to `fastp/`.

- [ ] **Step 4: Run tests and verify GREEN**

Run:

```powershell
& 'E:\Git\Git\bin\bash.exe' -n run_fastp_batch.sh
& 'E:\Git\Git\bin\bash.exe' tests/test_run_fastp_batch.sh
```

Expected: syntax exit 0 and test output ending in `PASS: run_fastp_batch.sh`.

- [ ] **Step 5: Commit implementation and tests**

Run:

```powershell
git add -- run_fastp_batch.sh tests/test_run_fastp_batch.sh
git commit -m "feat: add GEO fastp batch filter"
```

Expected: one commit containing only the script and its test.

### Task 3: Deploy without executing production filtering

**Files:**
- Copy: `run_fastp_batch.sh`
- Remote target: `/home/guanguiwen/10T4/murinematchtx_GEO/run_fastp_batch.sh`

**Interfaces:**
- Consumes: the locally verified script.
- Produces: one executable remote script whose SHA-256 matches the local source.

- [ ] **Step 1: Revalidate the SSH alias and remote target**

Run `ssh -G R540_2`, port reachability, the batch-mode SSH identity check, and:

```powershell
ssh -o BatchMode=yes R540_2 "test -d /home/guanguiwen/10T4/murinematchtx_GEO && test -f /home/guanguiwen/10T4/murinematchtx_GEO/srr_manifest.tsv"
```

Expected: exit 0.

- [ ] **Step 2: Copy exactly one file and set executable mode**

Run one scp command from the local script to the explicit target, then use `chmod 750` on that exact remote file. Do not invoke the script.

- [ ] **Step 3: Verify deployment**

Compare local and remote SHA-256 values, run remote `bash -n`, confirm the remote file mode, and run fastp 1.0.1 with 60 threads against empty standard input and reports directed to `/dev/null`.

Expected: hashes identical, syntax exit 0, mode executable, fastp smoke exit 0, and no `fastp/` output directories created under the production dataset root.

