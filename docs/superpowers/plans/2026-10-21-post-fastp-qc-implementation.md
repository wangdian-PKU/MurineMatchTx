# Post-fastp QC Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build, test, and deploy one safe Bash script that runs FastQC and MultiQC on all 109 fastp-filtered FASTQ files without starting the production QC run.

**Architecture:** The script validates the fixed production inventory, collects FASTQ files only from completed `fastp/` directories, and writes FastQC plus MultiQC output into a same-filesystem staging directory. It publishes `post_fastp_qc_results/` by atomic rename only after every expected report exists.

**Tech Stack:** Bash 5, FastQC 0.12.1, MultiQC 1.31, Git for Windows Bash for automated tests, PowerShell/scp/ssh for deployment.

## Global Constraints

- Local authoritative repository: `I:\R\Transcript similarity analysis\MurineMatchTx_remote_clean`.
- Remote project root: `/home/guanguiwen/10T4/murinematchtx_GEO` on `R540_2`.
- FastQC binary: `/home/guanguiwen/.conda/envs/RNAseq/bin/fastqc`.
- MultiQC binary: `/home/guanguiwen/.conda/envs/RNAseq/bin/multiqc`.
- Production inventory: exactly 66 `fastp/` directories with `.complete` and exactly 109 non-empty `*.fastq.gz` files.
- Input scope: only `<癌种>/<GSE>/<SRR>/fastp/*.fastq.gz`; never scan the raw `fastq/` directories.
- Final output: `post_fastp_qc_results/` with separate `fastqc/`, `multiqc/`, `fastqc.log`, and `multiqc.log` entries.
- Never overwrite, move, rename, or delete input data or existing QC results.
- Do not start the production FastQC or MultiQC run.

---

### Task 1: Implement the tested post-fastp QC script

**Files:**
- Create: `tests/test_run_post_fastp_qc.sh`
- Create: `run_post_fastp_qc.sh`

**Interfaces:**
- Consumes: `BASE_DIR`, `FASTQC_BIN`, and `MULTIQC_BIN` environment variables; production defaults point to the script directory and verified RNAseq-environment binaries.
- Produces: `<BASE_DIR>/post_fastp_qc_results/` only after complete FastQC and MultiQC validation; otherwise exits non-zero and retains the staging directory for diagnosis.

- [ ] **Step 1: Write the failing behavior test**

Create `tests/test_run_post_fastp_qc.sh` with this complete content:

```bash
#!/usr/bin/env bash
set -Eeuo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="${REPO_ROOT}/run_post_fastp_qc.sh"

if [[ ! -f "$SCRIPT" ]]; then
  printf 'FAIL: run_post_fastp_qc.sh is missing\n' >&2
  exit 1
fi

TEST_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/post-fastp-qc-test.XXXXXX")"
cleanup() {
  if [[ -n "${TEST_ROOT:-}" && -d "$TEST_ROOT" && "$TEST_ROOT" == *post-fastp-qc-test.* ]]; then
    find "$TEST_ROOT" -depth -delete
  fi
}
trap cleanup EXIT

DATA_ROOT="${TEST_ROOT}/data"
FAKE_FASTQC="${TEST_ROOT}/fake-fastqc"
FAKE_MULTIQC="${TEST_ROOT}/fake-multiqc"
FAIL_TOOL="${TEST_ROOT}/fail-tool"
CALL_LOG="${TEST_ROOT}/calls.log"
export CALL_LOG

mkdir -p "${DATA_ROOT}/1.COAD" "${DATA_ROOT}/2.LUAD" "${DATA_ROOT}/3.BRCA"
printf 'cancer\tgse\tsrr\tgsm\tlocal_source_column\tpublic_sample_title\n' \
  > "${DATA_ROOT}/srr_manifest.tsv"

for i in $(seq 1 66); do
  srr="$(printf 'SRR%08d' "$i")"
  srr_dir="${DATA_ROOT}/1.COAD/GSETEST/${srr}"
  mkdir -p "${srr_dir}/fastq" "${srr_dir}/fastp"
  printf 'raw-%s\n' "$srr" > "${srr_dir}/fastq/${srr}.fastq.gz"
  : > "${srr_dir}/fastp/.complete"

  if (( i <= 43 )); then
    printf 'filtered-r1-%s\n' "$srr" > "${srr_dir}/fastp/${srr}_1.fastq.gz"
    printf 'filtered-r2-%s\n' "$srr" > "${srr_dir}/fastp/${srr}_2.fastq.gz"
  else
    printf 'filtered-se-%s\n' "$srr" > "${srr_dir}/fastp/${srr}.fastq.gz"
  fi
done

cat > "$FAKE_FASTQC" <<'FAKE_FASTQC'
#!/usr/bin/env bash
set -Eeuo pipefail

threads=''
outdir=''
inputs=()
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --threads)
      threads="$2"
      shift 2
      ;;
    --outdir)
      outdir="$2"
      shift 2
      ;;
    *)
      inputs+=("$1")
      shift
      ;;
  esac
done

printf 'fastqc threads=%q outdir=%q' "$threads" "$outdir" >> "$CALL_LOG"
printf ' input=%q' "${inputs[@]}" >> "$CALL_LOG"
printf '\n' >> "$CALL_LOG"

[[ "$threads" == '80' ]]
[[ -d "$outdir" ]]
for input in "${inputs[@]}"; do
  name="${input##*/}"
  name="${name%.fastq.gz}"
  printf '<html>%s</html>\n' "$name" > "${outdir}/${name}_fastqc.html"
  printf 'zip-%s\n' "$name" > "${outdir}/${name}_fastqc.zip"
done
FAKE_FASTQC

cat > "$FAKE_MULTIQC" <<'FAKE_MULTIQC'
#!/usr/bin/env bash
set -Eeuo pipefail

outdir=''
filename=''
source_dir=''
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --outdir)
      outdir="$2"
      shift 2
      ;;
    --filename)
      filename="$2"
      shift 2
      ;;
    *)
      source_dir="$1"
      shift
      ;;
  esac
done

printf 'multiqc source=%q outdir=%q filename=%q\n' \
  "$source_dir" "$outdir" "$filename" >> "$CALL_LOG"
[[ -d "$source_dir" ]]
mkdir -p "${outdir}/multiqc_data"
printf '<html>multiqc</html>\n' > "${outdir}/${filename}"
FAKE_MULTIQC

cat > "$FAIL_TOOL" <<'FAIL_TOOL'
#!/usr/bin/env bash
exit 42
FAIL_TOOL

chmod 700 "$FAKE_FASTQC" "$FAKE_MULTIQC" "$FAIL_TOOL"

find "$DATA_ROOT" -type f -name '*.fastq.gz' -print0 \
  | sort -z \
  | xargs -0 sha256sum > "${TEST_ROOT}/before.sha256"

BASE_DIR="$DATA_ROOT" \
FASTQC_BIN="$FAKE_FASTQC" \
MULTIQC_BIN="$FAKE_MULTIQC" \
bash "$SCRIPT"

OUT="${DATA_ROOT}/post_fastp_qc_results"
[[ -d "${OUT}/fastqc" ]]
[[ -d "${OUT}/multiqc/multiqc_data" ]]
[[ -s "${OUT}/multiqc/multiqc_report.html" ]]
[[ -f "${OUT}/fastqc.log" ]]
[[ -f "${OUT}/multiqc.log" ]]
[[ "$(find "${OUT}/fastqc" -maxdepth 1 -type f -name '*_fastqc.html' | wc -l)" -eq 109 ]]
[[ "$(find "${OUT}/fastqc" -maxdepth 1 -type f -name '*_fastqc.zip' | wc -l)" -eq 109 ]]
[[ "$(wc -l < "$CALL_LOG")" -eq 2 ]]
grep -Fq 'fastqc threads=80' "$CALL_LOG"
grep -Fq 'multiqc source=' "$CALL_LOG"
! grep -Fq '/fastq/' "$CALL_LOG"

find "$DATA_ROOT" -type f -name '*.fastq.gz' -print0 \
  | sort -z \
  | xargs -0 sha256sum > "${TEST_ROOT}/after.sha256"
cmp "${TEST_ROOT}/before.sha256" "${TEST_ROOT}/after.sha256"

calls_before="$(wc -l < "$CALL_LOG")"
if BASE_DIR="$DATA_ROOT" FASTQC_BIN="$FAKE_FASTQC" MULTIQC_BIN="$FAKE_MULTIQC" \
  bash "$SCRIPT" >/dev/null 2>&1; then
  printf 'FAIL: existing output directory was accepted\n' >&2
  exit 1
fi
[[ "$(wc -l < "$CALL_LOG")" -eq "$calls_before" ]]

mv -- "$OUT" "${TEST_ROOT}/successful-output"
marker="$(find "$DATA_ROOT" -type f -path '*/fastp/.complete' | sort | head -n 1)"
rm -- "$marker"
if BASE_DIR="$DATA_ROOT" FASTQC_BIN="$FAKE_FASTQC" MULTIQC_BIN="$FAKE_MULTIQC" \
  bash "$SCRIPT" >/dev/null 2>&1; then
  printf 'FAIL: missing fastp completion marker was accepted\n' >&2
  exit 1
fi
[[ ! -e "$OUT" ]]
: > "$marker"

if BASE_DIR="$DATA_ROOT" FASTQC_BIN="$FAIL_TOOL" MULTIQC_BIN="$FAKE_MULTIQC" \
  bash "$SCRIPT" >/dev/null 2>&1; then
  printf 'FAIL: FastQC failure was accepted\n' >&2
  exit 1
fi
[[ ! -e "$OUT" ]]

if BASE_DIR="$DATA_ROOT" FASTQC_BIN="$FAKE_FASTQC" MULTIQC_BIN="$FAIL_TOOL" \
  bash "$SCRIPT" >/dev/null 2>&1; then
  printf 'FAIL: MultiQC failure was accepted\n' >&2
  exit 1
fi
[[ ! -e "$OUT" ]]

printf 'PASS: run_post_fastp_qc.sh\n'
```

- [ ] **Step 2: Run the test and verify RED**

Run:

```powershell
& 'E:\Git\Git\bin\bash.exe' tests/test_run_post_fastp_qc.sh
```

Expected: non-zero exit with `FAIL: run_post_fastp_qc.sh is missing`.

- [ ] **Step 3: Write the minimal production script**

Create `run_post_fastp_qc.sh` with this complete content:

```bash
#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${BASE_DIR:-$SCRIPT_DIR}"
FASTQC_BIN="${FASTQC_BIN:-/home/guanguiwen/.conda/envs/RNAseq/bin/fastqc}"
MULTIQC_BIN="${MULTIQC_BIN:-/home/guanguiwen/.conda/envs/RNAseq/bin/multiqc}"
THREADS=80
EXPECTED_FASTP_DIRS=66
EXPECTED_FASTQS=109
OUTPUT_DIR="${BASE_DIR}/post_fastp_qc_results"

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

[[ -d "$BASE_DIR" ]] || die "base directory does not exist: $BASE_DIR"
[[ -r "${BASE_DIR}/srr_manifest.tsv" ]] \
  || die "manifest is not readable: ${BASE_DIR}/srr_manifest.tsv"
[[ -x "$FASTQC_BIN" ]] || die "FastQC is not executable: $FASTQC_BIN"
[[ -x "$MULTIQC_BIN" ]] || die "MultiQC is not executable: $MULTIQC_BIN"
[[ ! -e "$OUTPUT_DIR" ]] || die "output already exists: $OUTPUT_DIR"

data_dirs=(
  "${BASE_DIR}/1.COAD"
  "${BASE_DIR}/2.LUAD"
  "${BASE_DIR}/3.BRCA"
)
for data_dir in "${data_dirs[@]}"; do
  [[ -d "$data_dir" ]] || die "data directory does not exist: $data_dir"
done

mapfile -d '' -t fastp_dirs < <(
  find "${data_dirs[@]}" \
    -mindepth 3 -maxdepth 3 -type d -name fastp -print0 \
    | sort -z
)
[[ "${#fastp_dirs[@]}" -eq "$EXPECTED_FASTP_DIRS" ]] \
  || die "expected $EXPECTED_FASTP_DIRS fastp directories, found ${#fastp_dirs[@]}"

for fastp_dir in "${fastp_dirs[@]}"; do
  [[ -f "${fastp_dir}/.complete" ]] \
    || die "fastp output is incomplete: $fastp_dir"
done

mapfile -d '' -t fastqs < <(
  find "${fastp_dirs[@]}" \
    -maxdepth 1 -type f -name '*.fastq.gz' -print0 \
    | sort -z
)
[[ "${#fastqs[@]}" -eq "$EXPECTED_FASTQS" ]] \
  || die "expected $EXPECTED_FASTQS filtered FASTQ files, found ${#fastqs[@]}"
for fastq in "${fastqs[@]}"; do
  [[ -s "$fastq" ]] || die "filtered FASTQ is empty: $fastq"
done

stage_dir="$(mktemp -d "${BASE_DIR}/.post_fastp_qc_results.XXXXXX")"
fastqc_dir="${stage_dir}/fastqc"
multiqc_dir="${stage_dir}/multiqc"
fastqc_log="${stage_dir}/fastqc.log"
multiqc_log="${stage_dir}/multiqc.log"
mkdir -p "$fastqc_dir" "$multiqc_dir"

printf 'Running FastQC on %s filtered FASTQ files with %s threads\n' \
  "${#fastqs[@]}" "$THREADS"
if ! "$FASTQC_BIN" \
  --threads "$THREADS" \
  --outdir "$fastqc_dir" \
  "${fastqs[@]}" \
  > "$fastqc_log" 2>&1; then
  die "FastQC failed; staged output retained at $stage_dir"
fi

fastqc_html_count="$(
  find "$fastqc_dir" -maxdepth 1 -type f -name '*_fastqc.html' | wc -l
)"
fastqc_zip_count="$(
  find "$fastqc_dir" -maxdepth 1 -type f -name '*_fastqc.zip' | wc -l
)"
[[ "$fastqc_html_count" -eq "$EXPECTED_FASTQS" ]] \
  || die "expected $EXPECTED_FASTQS FastQC HTML reports, found $fastqc_html_count; staged output retained at $stage_dir"
[[ "$fastqc_zip_count" -eq "$EXPECTED_FASTQS" ]] \
  || die "expected $EXPECTED_FASTQS FastQC ZIP reports, found $fastqc_zip_count; staged output retained at $stage_dir"

printf 'Running MultiQC\n'
if ! "$MULTIQC_BIN" \
  "$fastqc_dir" \
  --outdir "$multiqc_dir" \
  --filename multiqc_report.html \
  > "$multiqc_log" 2>&1; then
  die "MultiQC failed; staged output retained at $stage_dir"
fi

[[ -s "${multiqc_dir}/multiqc_report.html" ]] \
  || die "MultiQC report is missing or empty; staged output retained at $stage_dir"
[[ -d "${multiqc_dir}/multiqc_data" ]] \
  || die "MultiQC data directory is missing; staged output retained at $stage_dir"

mv -- "$stage_dir" "$OUTPUT_DIR"
printf 'Post-fastp QC finished: %s\n' "$OUTPUT_DIR"
```

- [ ] **Step 4: Run syntax and behavior tests and verify GREEN**

Run:

```powershell
& 'E:\Git\Git\bin\bash.exe' -n run_post_fastp_qc.sh
& 'E:\Git\Git\bin\bash.exe' -n tests/test_run_post_fastp_qc.sh
& 'E:\Git\Git\bin\bash.exe' tests/test_run_post_fastp_qc.sh
```

Expected: both syntax checks exit 0; the behavior test ends with `PASS: run_post_fastp_qc.sh`.

- [ ] **Step 5: Review the exact change and commit**

Run:

```powershell
git update-index --add --chmod=+x run_post_fastp_qc.sh tests/test_run_post_fastp_qc.sh
git diff --cached --check
git diff --cached -- run_post_fastp_qc.sh tests/test_run_post_fastp_qc.sh
git commit -m "feat: add post-fastp quality control"
```

Expected: one commit containing only the production script and its behavior test.

---

### Task 2: Deploy the script without executing QC

**Files:**
- Copy: `run_post_fastp_qc.sh`
- Remote target: `/home/guanguiwen/10T4/murinematchtx_GEO/run_post_fastp_qc.sh`

**Interfaces:**
- Consumes: the locally tested and committed script from Task 1.
- Produces: one executable remote script whose SHA-256 matches the local source; does not produce QC results.

- [ ] **Step 1: Revalidate SSH and the remote production inventory**

Run:

```powershell
$OutputEncoding = [Console]::OutputEncoding = [System.Text.UTF8Encoding]::new()
ssh -G R540_2 | Select-String '^(hostname|port|user|identityfile) '
Test-NetConnection frp6.ccszxc.xin -Port 42246 -InformationLevel Quiet
ssh -o BatchMode=yes -o ConnectTimeout=10 R540_2 "hostname; whoami; pwd"
ssh -o BatchMode=yes R540_2 "test -d /home/guanguiwen/10T4/murinematchtx_GEO && test -f /home/guanguiwen/10T4/murinematchtx_GEO/srr_manifest.tsv"
ssh -o BatchMode=yes R540_2 "find /home/guanguiwen/10T4/murinematchtx_GEO -mindepth 5 -maxdepth 5 -type f -path '*/fastp/.complete' | wc -l"
ssh -o BatchMode=yes R540_2 "find /home/guanguiwen/10T4/murinematchtx_GEO -mindepth 5 -maxdepth 5 -type f -path '*/fastp/*.fastq.gz' | wc -l"
```

Expected: alias resolves to `frp6.ccszxc.xin:42246` as `guanguiwen`; connection succeeds; the final two counts are `66` and `109`.

- [ ] **Step 2: Confirm the exact remote target is unused, then copy one file**

Run:

```powershell
ssh -o BatchMode=yes R540_2 "test ! -e /home/guanguiwen/10T4/murinematchtx_GEO/run_post_fastp_qc.sh"
scp 'I:\R\Transcript similarity analysis\MurineMatchTx_remote_clean\run_post_fastp_qc.sh' 'R540_2:/home/guanguiwen/10T4/murinematchtx_GEO/run_post_fastp_qc.sh'
ssh -o BatchMode=yes R540_2 "chmod 750 /home/guanguiwen/10T4/murinematchtx_GEO/run_post_fastp_qc.sh"
```

Expected: the non-existence check and all transfer commands exit 0. If the target already exists, stop rather than overwrite it.

- [ ] **Step 3: Verify deployment without running production QC**

Run:

```powershell
(Get-FileHash -LiteralPath 'I:\R\Transcript similarity analysis\MurineMatchTx_remote_clean\run_post_fastp_qc.sh' -Algorithm SHA256).Hash.ToLowerInvariant()
ssh -o BatchMode=yes R540_2 "sha256sum /home/guanguiwen/10T4/murinematchtx_GEO/run_post_fastp_qc.sh"
ssh -o BatchMode=yes R540_2 "bash -n /home/guanguiwen/10T4/murinematchtx_GEO/run_post_fastp_qc.sh"
ssh -o BatchMode=yes R540_2 "stat -c %a /home/guanguiwen/10T4/murinematchtx_GEO/run_post_fastp_qc.sh"
ssh -o BatchMode=yes R540_2 "test ! -e /home/guanguiwen/10T4/murinematchtx_GEO/post_fastp_qc_results"
```

Expected: SHA-256 values are identical, syntax exits 0, mode is `750`, and the final check confirms that production QC has not started.
