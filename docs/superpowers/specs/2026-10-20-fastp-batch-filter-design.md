# R540_2 GEO FASTQ 批量 fastp 过滤设计

## 目标与边界

为 `/home/guanguiwen/10T4/murinematchtx_GEO` 下 3 个癌种、9 个 GSE、66 个 SRR 生成一个批量过滤脚本。脚本只读取各 SRR 的 `fastq/` 原始数据，不覆盖、不移动、不删除原始 FASTQ，也不自动开始全量过滤。

fastp 使用已验证的绝对路径：

`/home/guanguiwen/.conda/envs/fastp/bin/fastp`（fastp 1.0.1）

## 数据集规则

| 癌种 | GSE | 布局 | 接头处理 | 特殊处理 |
| --- | --- | --- | --- | --- |
| COAD | GSE245724 | PE | 关闭 | 无 |
| COAD | GSE263502 | SE | 关闭 | 无 |
| COAD | GSE285835 | PE | 启用 PE 自动检测 | 无 |
| LUAD | GSE109629 | PE | 关闭 | 无 |
| LUAD | GSE163087 | SE | 关闭 | 无 |
| LUAD | GSE70046 | PE | 关闭 | R1、R2 各精确去除 3′ 末端 1 bp |
| BRCA | GSE118802 | PE | 启用 PE 自动检测 | 无 |
| BRCA | GSE159203 | PE | 关闭 | 无 |
| BRCA | GSE262136 | PE | 关闭 | 无 |

GSE118802 的全部 12 个 FASTQ 在 FastQC 中出现 Adapter Content FAIL；GSE285835 的 SRR31881982 两端出现 Adapter Content WARN。GSE70046 的全部 10 个 FASTQ 出现末端质量/N 含量 FAIL，采用固定末端 1 bp 修剪，不对其他数据集启用通用滑窗修剪。

## 固定 fastp 参数

- `--thread 60`
- `--qualified_quality_phred 20`
- `--unqualified_percent_limit 40`
- `--n_base_limit 5`
- `--length_required 30`
- `--disable_trim_poly_g`
- 不启用 `cut_front`、`cut_tail`、`cut_right`，保留 5′ 前 10 bp。
- 不启用去重复、低复杂度过滤、polyX 修剪、碱基纠错或 read 合并。
- PE 数据任一端未通过过滤时按 fastp 默认行为丢弃整对，不生成额外 unpaired 文件。

接头阳性的 PE 数据使用 `--detect_adapter_for_pe`；其余数据显式使用 `--disable_adapter_trimming`。

## 输入、输出与安全行为

输入结构：

`<癌种>/<GSE>/<SRR>/fastq/*.fastq.gz`

输出结构：

`<癌种>/<GSE>/<SRR>/fastp/`

输出 FASTQ 保留原始文件名，并在独立目录内生成 `fastp.json`、`fastp.html`、`fastp.log` 和 `.complete`。每个 SRR 先写入同一 SRR 目录下的临时目录，fastp 成功后再原子移动为 `fastp/`。

- `fastp/.complete` 已存在：跳过该 SRR。
- `fastp/` 已存在但无 `.complete`：立即报错，不覆盖、不自动删除。
- 输入缺失、单/双端文件数量不符合规则、报告或输出为空：立即报错并停止。
- 全部样本顺序执行；任一时刻只有一个 fastp 进程使用 60 线程。

## 验证标准

实现前先建立最小测试夹具并验证测试因脚本缺失而失败；实现后必须满足：

1. SE、普通 PE、接头阳性 PE、GSE70046 四类参数分派正确。
2. 所有命令都包含 Q20、`-u 40`、N≤5、最短 30 bp 和 60 线程。
3. 所有命令均不包含任何 5′ 修剪、去重复或低复杂度过滤参数。
4. 原始输入哈希在测试前后不变；已有完整输出被跳过，已有不完整输出导致非零退出。
5. Shell 语法检查、自动化测试和小型真实 fastp 冒烟测试全部退出码为 0。

交付脚本只负责批量过滤，不自动执行 189 GiB 原始数据的全量任务。
