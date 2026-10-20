#!/usr/bin/env bash
set -Eeuo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPT="${REPO_ROOT}/run_fastp_batch.sh"

if [[ ! -f "$SCRIPT" ]]; then
  printf 'FAIL: run_fastp_batch.sh is missing\n' >&2
  exit 1
fi

TEST_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/fastp-batch-test.XXXXXX")"
cleanup() {
  if [[ -n "${TEST_ROOT:-}" && -d "$TEST_ROOT" && "$TEST_ROOT" == *fastp-batch-test.* ]]; then
    find "$TEST_ROOT" -depth -delete
  fi
}
trap cleanup EXIT

DATA_ROOT="${TEST_ROOT}/data"
MANIFEST="${TEST_ROOT}/manifest.tsv"
FAKE_FASTP="${TEST_ROOT}/fake-fastp"
CALL_LOG="${TEST_ROOT}/calls.log"
export CALL_LOG

mkdir -p \
  "${DATA_ROOT}/1.COAD/GSE263502/SRR28591489/fastq" \
  "${DATA_ROOT}/1.COAD/GSE245724/SRR26425828/fastq" \
  "${DATA_ROOT}/1.COAD/GSE285835/SRR31881982/fastq" \
  "${DATA_ROOT}/2.LUAD/GSE163087/SRR13246552/fastq" \
  "${DATA_ROOT}/2.LUAD/GSE70046/SRR2071501/fastq"

printf 'se-original\n' > "${DATA_ROOT}/1.COAD/GSE263502/SRR28591489/fastq/SRR28591489.fastq.gz"
printf 'pe-normal-r1\n' > "${DATA_ROOT}/1.COAD/GSE245724/SRR26425828/fastq/SRR26425828_1.fastq.gz"
printf 'pe-normal-r2\n' > "${DATA_ROOT}/1.COAD/GSE245724/SRR26425828/fastq/SRR26425828_2.fastq.gz"
printf 'pe-adapter-r1\n' > "${DATA_ROOT}/1.COAD/GSE285835/SRR31881982/fastq/SRR31881982_1.fastq.gz"
printf 'pe-adapter-r2\n' > "${DATA_ROOT}/1.COAD/GSE285835/SRR31881982/fastq/SRR31881982_2.fastq.gz"
printf 'se-read1-original\n' > "${DATA_ROOT}/2.LUAD/GSE163087/SRR13246552/fastq/SRR13246552_1.fastq.gz"
printf 'pe-tail-r1\n' > "${DATA_ROOT}/2.LUAD/GSE70046/SRR2071501/fastq/SRR2071501_1.fastq.gz"
printf 'pe-tail-r2\n' > "${DATA_ROOT}/2.LUAD/GSE70046/SRR2071501/fastq/SRR2071501_2.fastq.gz"

printf '%s\n' \
  $'cancer\tgse\tsrr\tgsm\tlocal_source_column\tpublic_sample_title' \
  $'1.COAD\tGSE263502\tSRR28591489\tGSM1\tSE\tSE sample' \
  $'1.COAD\tGSE245724\tSRR26425828\tGSM2\tPE\tPE sample' \
  $'1.COAD\tGSE285835\tSRR31881982\tGSM3\tADAPTER\tAdapter sample' \
  $'2.LUAD\tGSE163087\tSRR13246552\tGSM5\tSE_R1\tSE _1 sample' \
  $'2.LUAD\tGSE70046\tSRR2071501\tGSM4\tTAIL\tTail sample' \
  > "$MANIFEST"

cat > "$FAKE_FASTP" <<'FAKE'
#!/usr/bin/env bash
set -Eeuo pipefail

printf '%q ' "$@" >> "$CALL_LOG"
printf '\n' >> "$CALL_LOG"

out1=''
out2=''
json=''
html=''
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --out1)
      out1="$2"
      shift 2
      ;;
    --out2)
      out2="$2"
      shift 2
      ;;
    --json)
      json="$2"
      shift 2
      ;;
    --html)
      html="$2"
      shift 2
      ;;
    *)
      shift
      ;;
  esac
done

printf 'filtered-r1\n' > "$out1"
if [[ -n "$out2" ]]; then
  printf 'filtered-r2\n' > "$out2"
fi
printf '{}\n' > "$json"
printf '<html></html>\n' > "$html"
FAKE
chmod 700 "$FAKE_FASTP"

find "$DATA_ROOT" -type f -name '*.fastq.gz' -print0 \
  | sort -z \
  | xargs -0 sha256sum > "${TEST_ROOT}/before.sha256"

BASE_DIR="$DATA_ROOT" FASTP_BIN="$FAKE_FASTP" bash "$SCRIPT" "$MANIFEST"

for complete in \
  "${DATA_ROOT}/1.COAD/GSE263502/SRR28591489/fastp/.complete" \
  "${DATA_ROOT}/1.COAD/GSE245724/SRR26425828/fastp/.complete" \
  "${DATA_ROOT}/1.COAD/GSE285835/SRR31881982/fastp/.complete" \
  "${DATA_ROOT}/2.LUAD/GSE163087/SRR13246552/fastp/.complete" \
  "${DATA_ROOT}/2.LUAD/GSE70046/SRR2071501/fastp/.complete"; do
  [[ -f "$complete" ]] || { printf 'FAIL: missing %s\n' "$complete" >&2; exit 1; }
done

[[ "$(wc -l < "$CALL_LOG")" -eq 5 ]]
while IFS= read -r call; do
  grep -Fq -- '--thread 60' <<< "$call"
  grep -Fq -- '--qualified_quality_phred 20' <<< "$call"
  grep -Fq -- '--unqualified_percent_limit 40' <<< "$call"
  grep -Fq -- '--n_base_limit 5' <<< "$call"
  grep -Fq -- '--length_required 30' <<< "$call"
  grep -Fq -- '--disable_trim_poly_g' <<< "$call"
  ! grep -Eq -- '--trim_front|--cut_front|--cut_tail|--cut_right|--dedup|--low_complexity_filter|--trim_poly_x|--correction|--merge' <<< "$call"
done < "$CALL_LOG"

[[ "$(grep -Fc -- '--detect_adapter_for_pe' "$CALL_LOG")" -eq 1 ]]
[[ "$(grep -Fc -- '--disable_adapter_trimming' "$CALL_LOG")" -eq 4 ]]
[[ "$(grep -Fc -- '--trim_tail1 1 --trim_tail2 1' "$CALL_LOG")" -eq 1 ]]
[[ "$(grep -Fc -- '--in2' "$CALL_LOG")" -eq 3 ]]

find "$DATA_ROOT" -type f -path '*/fastq/*.fastq.gz' -print0 \
  | sort -z \
  | xargs -0 sha256sum > "${TEST_ROOT}/after.sha256"
cmp "${TEST_ROOT}/before.sha256" "${TEST_ROOT}/after.sha256"

calls_before="$(wc -l < "$CALL_LOG")"
BASE_DIR="$DATA_ROOT" FASTP_BIN="$FAKE_FASTP" bash "$SCRIPT" "$MANIFEST"
[[ "$(wc -l < "$CALL_LOG")" -eq "$calls_before" ]]

rm -- "${DATA_ROOT}/2.LUAD/GSE70046/SRR2071501/fastp/.complete"
if BASE_DIR="$DATA_ROOT" FASTP_BIN="$FAKE_FASTP" bash "$SCRIPT" "$MANIFEST" >/dev/null 2>&1; then
  printf 'FAIL: incomplete fastp directory was accepted\n' >&2
  exit 1
fi

printf 'PASS: run_fastp_batch.sh\n'
