# fastp 过滤后批量 QC 设计

## 目标与边界

为 `/home/guanguiwen/10T4/murinematchtx_GEO` 下已经完成 fastp 过滤的 66 个 SRR 提供一个批量 QC 脚本。脚本只读取各 SRR 的 `fastp/` 过滤后 FASTQ，不读取原始 `fastq/`，不修改任何输入，也不自动开始正式 QC。

脚本放在项目根目录，命名为 `run_post_fastp_qc.sh`。运行环境使用已经验证的工具：

- `/home/guanguiwen/.conda/envs/RNAseq/bin/fastqc`（FastQC 0.12.1）
- `/home/guanguiwen/.conda/envs/RNAseq/bin/multiqc`（MultiQC 1.31）

## 输入与完整性检查

输入范围严格限定为：

`<癌种>/<GSE>/<SRR>/fastp/*.fastq.gz`

脚本运行前必须确认：

1. 项目根目录及 `srr_manifest.tsv` 存在。
2. 66 个 SRR 均存在 `fastp/.complete`。
3. 共发现 109 个非空的过滤后 FASTQ。
4. FastQC 和 MultiQC 的绝对路径均可执行。
5. 正式输出目录尚不存在。

任一检查失败时立即退出，不创建正式结果目录。固定数量检查用于防止 fastp 未完成、文件缺失或误把原始 FASTQ 纳入 QC。

## 处理流程与输出

脚本在项目根目录下创建同文件系统的临时目录，按稳定排序收集 109 个过滤后 FASTQ，并以 80 个线程运行 FastQC。随后 MultiQC 只扫描本轮 FastQC 输出，不扫描项目其他目录。

正式输出目录为：

```text
post_fastp_qc_results/
├── fastqc/
│   ├── *_fastqc.html
│   └── *_fastqc.zip
├── multiqc/
│   ├── multiqc_report.html
│   └── multiqc_data/
├── fastqc.log
└── multiqc.log
```

FastQC 和 MultiQC 都成功，且生成 109 份 FastQC HTML、109 份 FastQC ZIP 和非空的 `multiqc_report.html` 后，脚本才将临时目录原子改名为 `post_fastp_qc_results/`。若执行失败，保留临时目录和日志以便诊断；不生成看似完整的正式目录。

## 安全与重复执行

- 输入 FASTQ 和既有 `fastqc_results/` 均保持只读。
- 不覆盖、不删除、不移动任何已有文件。
- `post_fastp_qc_results/` 已存在时立即报错，由使用者决定如何处理旧结果。
- 使用 Bash 严格模式；路径全部加引号；通过空字符分隔文件名。
- 用户自行执行脚本，本次交付不启动正式 QC。

## 验证标准

实现前先建立最小测试夹具，并确认测试因脚本缺失而失败。实现后必须满足：

1. 只把 `fastp/` 内的过滤后 FASTQ 传给 FastQC，原始 `fastq/` 不会被扫描。
2. 缺少 `.complete`、FASTQ 为空、数量不符或输出目录已存在时以非零状态退出。
3. FastQC 和 MultiQC 使用预期参数，MultiQC 只汇总本轮 FastQC 目录。
4. 工具失败时不创建正式输出目录，成功时才原子发布完整结果。
5. Shell 语法检查和自动化行为测试均退出码为 0。
