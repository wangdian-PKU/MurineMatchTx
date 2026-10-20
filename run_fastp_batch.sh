#!/usr/bin/env bash
set -Eeuo pipefail

shopt -s nullglob

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${BASE_DIR:-$SCRIPT_DIR}"
FASTP_BIN="${FASTP_BIN:-/home/guanguiwen/.conda/envs/fastp/bin/fastp}"
MANIFEST="${1:-${SCRIPT_DIR}/srr_manifest.tsv}"
THREADS=60
CANONICAL_MANIFEST="${SCRIPT_DIR}/srr_manifest.tsv"

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

[[ -x "$FASTP_BIN" ]] || die "fastp is not executable: $FASTP_BIN"
[[ -r "$MANIFEST" ]] || die "manifest is not readable: $MANIFEST"
[[ -d "$BASE_DIR" ]] || die "base directory does not exist: $BASE_DIR"

row_count="$(awk -F '\t' 'NR > 1 && NF { n++ } END { print n + 0 }' "$MANIFEST")"
if [[ "$MANIFEST" == "$CANONICAL_MANIFEST" && "$row_count" -ne 66 ]]; then
  die "expected 66 data rows in $MANIFEST, found $row_count"
fi
[[ "$row_count" -gt 0 ]] || die "manifest contains no data rows: $MANIFEST"

declare -A seen_srr=()
processed=0
skipped=0

while IFS=$'\t' read -r cancer gse srr gsm local_source_column public_sample_title; do
  if [[ "$cancer" == "cancer" || -z "$cancer" ]]; then
    continue
  fi

  public_sample_title="${public_sample_title%$'\r'}"
  [[ "$cancer" =~ ^[123]\.(COAD|LUAD|BRCA)$ ]] \
    || die "invalid cancer value in manifest: $cancer"
  [[ "$gse" =~ ^GSE[0-9]+$ ]] || die "invalid GSE value in manifest: $gse"
  [[ "$srr" =~ ^SRR[0-9]+$ ]] || die "invalid SRR value in manifest: $srr"
  [[ -z "${seen_srr[$srr]:-}" ]] || die "duplicate SRR in manifest: $srr"
  seen_srr[$srr]=1

  layout=''
  adapter='disable'
  tail_trim=0
  case "${cancer}/${gse}" in
    1.COAD/GSE263502|2.LUAD/GSE163087)
      layout='SE'
      ;;
    1.COAD/GSE285835|3.BRCA/GSE118802)
      layout='PE'
      adapter='detect'
      ;;
    2.LUAD/GSE70046)
      layout='PE'
      tail_trim=1
      ;;
    1.COAD/GSE245724|2.LUAD/GSE109629|3.BRCA/GSE159203|3.BRCA/GSE262136)
      layout='PE'
      ;;
    *)
      die "unsupported dataset in manifest: ${cancer}/${gse}"
      ;;
  esac

  srr_dir="${BASE_DIR}/${cancer}/${gse}/${srr}"
  fastq_dir="${srr_dir}/fastq"
  output_dir="${srr_dir}/fastp"
  single_in="${fastq_dir}/${srr}.fastq.gz"
  read1_in="${fastq_dir}/${srr}_1.fastq.gz"
  read2_in="${fastq_dir}/${srr}_2.fastq.gz"

  [[ -d "$fastq_dir" ]] || die "FASTQ directory does not exist: $fastq_dir"
  if [[ -f "${output_dir}/.complete" ]]; then
    printf '[%s/%s] already complete: %s\n' "$cancer" "$gse" "$srr"
    skipped=$((skipped + 1))
    continue
  fi
  [[ ! -e "$output_dir" ]] \
    || die "incomplete output directory exists; inspect it before rerunning: $output_dir"

  input1=''
  input2=''
  if [[ "$layout" == 'SE' ]]; then
    if [[ -s "$single_in" && ! -e "$read1_in" && ! -e "$read2_in" ]]; then
      input1="$single_in"
    elif [[ -s "$read1_in" && ! -e "$single_in" && ! -e "$read2_in" ]]; then
      input1="$read1_in"
    else
      die "expected exactly one SE FASTQ (${srr}.fastq.gz or ${srr}_1.fastq.gz): $fastq_dir"
    fi
  else
    [[ -s "$read1_in" && -s "$read2_in" && ! -e "$single_in" ]] \
      || die "expected one non-empty PE FASTQ pair for $srr: $fastq_dir"
    input1="$read1_in"
    input2="$read2_in"
  fi

  stage_dir="$(mktemp -d "${srr_dir}/.fastp.XXXXXX")"
  output1="${stage_dir}/$(basename -- "$input1")"
  output2=''
  json_report="${stage_dir}/fastp.json"
  html_report="${stage_dir}/fastp.html"
  log_file="${stage_dir}/fastp.log"

  args=(
    --in1 "$input1"
    --out1 "$output1"
    --json "$json_report"
    --html "$html_report"
    --report_title "${cancer} ${gse} ${srr}"
    --thread "$THREADS"
    --qualified_quality_phred 20
    --unqualified_percent_limit 40
    --n_base_limit 5
    --length_required 30
    --disable_trim_poly_g
    --dont_overwrite
  )

  if [[ "$layout" == 'PE' ]]; then
    output2="${stage_dir}/$(basename -- "$input2")"
    args+=(--in2 "$input2" --out2 "$output2")
  fi

  if [[ "$adapter" == 'detect' ]]; then
    args+=(--detect_adapter_for_pe)
  else
    args+=(--disable_adapter_trimming)
  fi

  if [[ "$tail_trim" -eq 1 ]]; then
    args+=(--trim_tail1 1 --trim_tail2 1)
  fi

  printf '[%s/%s] fastp: %s (%s)\n' "$cancer" "$gse" "$srr" "$layout"
  if ! "$FASTP_BIN" "${args[@]}" > "$log_file" 2>&1; then
    die "fastp failed for $srr; staged output retained at $stage_dir"
  fi

  [[ -s "$output1" ]] || die "fastp produced no R1/SE output for $srr: $output1"
  if [[ "$layout" == 'PE' ]]; then
    [[ -s "$output2" ]] || die "fastp produced no R2 output for $srr: $output2"
  fi
  [[ -s "$json_report" ]] || die "fastp JSON report is empty for $srr"
  [[ -s "$html_report" ]] || die "fastp HTML report is empty for $srr"

  touch -- "${stage_dir}/.complete"
  mv -- "$stage_dir" "$output_dir"
  processed=$((processed + 1))
done < "$MANIFEST"

if [[ "${#seen_srr[@]}" -ne "$row_count" ]]; then
  die "processed manifest row count mismatch: expected $row_count, saw ${#seen_srr[@]}"
fi

printf 'fastp batch finished: processed=%s skipped=%s total=%s\n' \
  "$processed" "$skipped" "$row_count"
