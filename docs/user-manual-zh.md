# phyluce Rust CLI 中文用户手册

本文档面向已经熟悉原版 `phyluce` 的用户，说明 Rust 版本如何安装、配置和使用。Rust 版本的目标是模仿原版脚本的工作方式，同时提供一个更统一的命令入口。

当前 Rust CLI 处于早期预览阶段。多数命令已经按原版脚本行为实现，但仍有少数命令依赖外部程序，或与原版存在有意的兼容差异。正式分析前，建议先用小数据或现有 fixture 运行一遍完整流程。

## 0. 阅读指南与适用范围

本手册只描述 `phyluce` Rust CLI。原版 PHYLUCE 的
[官方文档](https://phyluce.readthedocs.io/en/latest/index.html)适合了解 UCE
方法背景和经典工作流，但其中的 Python/Conda 安装命令不能直接用于本项目。
Rust 版本的实际参数以 `phyluce --help` 和各级子命令帮助为准。

Rust CLI 主要覆盖以下环节：

- 调用外部 assembler，并整理各样本 contigs。
- 将 contigs 与 UCE probes 匹配，过滤重复或疑似旁系同源命中。
- 生成完整或不完整 locus 集合并提取序列。
- 使用 MAFFT 比对，随后修剪、过滤、统计和拼接 alignment。
- 处理 probe/genome、gene-tree、NCBI 提交和常用格式转换任务。

以下工作不属于本 CLI 的完整职责：

- 原始 reads 的接头去除、低质量碱基清理和实验批次质量评估。
- MAFFT、LASTZ、SPAdes、Gblocks、trimAl、RAxML、Snakemake 等外部程序
  本身的安装和算法参数选择。
- 使用 IQ-TREE、RAxML、ASTRAL 等完成最终系统发育推断。

推荐阅读路线：

| 目标 | 阅读章节 |
| --- | --- |
| 首次安装并跑通流程 | 2 -> 3 -> 4 -> 6 -> 7 -> 17 |
| 从原版 `phyluce_*` 脚本迁移 | 1 -> 5 -> 14 -> 16 |
| 复现实验或报告错误 | 4 -> 15 -> 18 -> 19 -> 20 |
| 只查询某个命令 | 6-13，并同时查看对应的 `--help` |

## 1. 命令形式

Rust 版本提供一个统一入口：

```bash
phyluce <domain> <command> [options]
```

例如，原版：

```bash
phyluce_align_convert_degen_bases --alignments in --output out
```

Rust 版本：

```bash
phyluce align convert-degen-bases --alignments in --output out
```

命令按功能分组：

| 分组 | 用途 |
| --- | --- |
| `config` | 查看和解析 `phyluce.conf` 配置 |
| `io` | FASTA/FASTQ 输入检查 |
| `external` | 检查外部程序路径 |
| `assembly` | assembly、contig 匹配、match-count 和 FASTA 提取 |
| `align` | alignment、trim、转换、过滤、拼接和统计 |
| `probe` | 探针设计、lastz 结果解析、BED/2bit 操作 |
| `utilities` | 原版 utilities 小工具 |
| `ncbi` | NCBI 提交前 FASTA 准备 |
| `genetrees` | gene tree 重命名和统计 |
| `workflow` | 调用 Snakemake workflow |

查看帮助：

```bash
phyluce --help
phyluce assembly --help
phyluce assembly match-contigs-to-probes --help
```

## 2. 构建与安装

### 2.1 环境要求

- 支持 Rust 2021 edition 的稳定版 Rust 工具链和 Cargo。
- 推荐 Linux 或 macOS 等类 Unix 环境。实际可移植性还取决于所调用的外部程序。
- 足够保存 reads、contigs、逐 locus alignment、LASTZ 中间结果和最终矩阵的
  磁盘空间。大型 LASTZ 与拼接任务还需要可用的临时目录空间。
- 只安装当前分析会调用的外部程序。Conda 可用于管理这些程序，但不是 Rust
  CLI 的安装方式。

### 2.2 构建

在仓库根目录运行：

```bash
cargo build -p phyluce-cli --release
```

构建后的程序位于：

```bash
target/release/phyluce
```

可以直接运行：

```bash
target/release/phyluce --help
```

也可以把它加入 `PATH`：

```bash
export PATH="/path/to/phyluce_rust/target/release:$PATH"
```

开发检查：

```bash
cargo check -p phyluce-cli
cargo test -p phyluce-io -p phyluce-assembly -p phyluce-cli
cargo clippy -p phyluce-cli --all-targets -- -D warnings
```

### 2.3 推荐项目目录

每个分析阶段使用独立输出目录，便于回溯、重跑和比较参数：

```text
project/
  raw-fastq/                 # 原始数据，只读保存
  clean-fastq/               # 接头和质量清理后的 reads
  assembly/
    contigs/
  probes/
    uce-probes.fasta
  uce-search/                # LASTZ 与 probe.matches.sqlite
  taxon-sets/                # taxon 配置和 complete/incomplete locus 列表
  alignments/
    untrimmed/
    trimmed/
    filtered/
  matrices/                  # NEXUS/PHYLIP 与 charsets
  stats/
  logs/
```

不要把输出写回输入目录，也不要在确认新结果前覆盖旧阶段。分析开始时记录
`phyluce --version`、外部程序版本和实际使用的配置文件；第 19 节给出完整的
复现清单。

## 3. 配置外部程序

Rust 版本沿用原版 `phyluce.conf` 的思路。配置中通常包含：

```ini
[binaries]
lastz:$CONDA/bin/lastz
mafft:$CONDA/bin/mafft
gblocks:$CONDA/bin/Gblocks
trimal:$CONDA/bin/trimal
raxml-ng:$CONDA/bin/raxml-ng
spades:$CONDA/bin/spades.py
velveth:$CONDA/bin/velveth
velvetg:$CONDA/bin/velvetg
abyss:$CONDA/bin/ABYSS
abyss-pe:$CONDA/bin/abyss-pe
samtools:$CONDA/bin/samtools
bcftools:$CONDA/bin/bcftools
snakemake:$CONDA/bin/Snakemake

[workflows]
mapping:$WORKFLOWS/mapping/Snakefile
correction:$WORKFLOWS/contig-correction/Snakefile
phasing:$WORKFLOWS/phasing/Snakefile
```

Rust 版本内嵌了上述默认配置，并在存在 `config/phyluce.conf` 时优先读取磁盘
版本，随后合并 `~/.phyluce.conf`。因此独立安装不需要保留源码目录。程序会
展开 `$CONDA`、`$WORKFLOWS`；也可通过 `PHYLUCE_CONFIG` 指定配置文件。

查看当前解析到的配置：

```bash
phyluce config inspect
```

查看某个外部程序路径：

```bash
phyluce config which --program binaries --binary lastz
phyluce config which --program binaries --binary mafft
```

检查外部程序是否可运行：

```bash
phyluce external check --program binaries --binary lastz
```

## 4. 日志

所有命令都支持全局日志参数：

```bash
phyluce --log-path log --verbosity INFO assembly get-fasta-lengths --input contigs.fasta
```

参数说明：

| 参数 | 说明 |
| --- | --- |
| `--log-path DIR` | 写入 `DIR/phyluce.log` |
| `--verbosity INFO|WARN|CRITICAL` | 控制日志等级，默认 `INFO` |

如果通过旧脚本名兼容方式运行，日志文件名会使用旧脚本名，例如：

```text
phyluce_align_convert_degen_bases.log
```

默认不写日志，以免改变 stdout/stderr 行为。

## 5. 原版脚本名兼容

Rust 二进制可识别全部 74 个原版可执行脚本名：如果把 `phyluce` 复制或软链接为旧脚本名，它会自动映射到新的分组命令。脚本名可映射不代表所有旧选项和中间文件格式完全等价，具体边界见“已知差异”。

例如：

```bash
ln -s target/release/phyluce phyluce_align_convert_degen_bases
./phyluce_align_convert_degen_bases --alignments in --output out
```

等价于：

```bash
phyluce align convert-degen-bases --alignments in --output out
```

例如：

```text
phyluce_align_convert_degen_bases
phyluce_assembly_get_fasta_lengths
phyluce_probe_easy_lastz
phyluce_genetrees_get_tree_counts
phyluce_utilities_combine_reads
phyluce_ncbi_chunk_fasta_for_ncbi
phyluce_workflow
```

所有原版入口都可以使用相同方式创建链接；新分组式 CLI 仍是推荐入口。

## 6. UCE 主流程

典型 UCE phylogenomics 流程如下：

```text
raw reads
  -> adapter/quality cleaning（外部工具）
  -> assembly
  -> contigs/*.contigs.fasta
  -> contigs 与 probes 做 LASTZ 匹配
  -> probe.matches.sqlite
  -> match-count config
  -> monolithic UCE FASTA
  -> 按 locus 比对
  -> trimming / filtering / concatenation / summary
```

进入流程前先确认：

1. reads 已完成接头去除和质量清理，并保留相应质控报告。
2. 样本名在 reads 目录、assembly 配置和 taxon 配置中完全一致；建议只使用
   字母、数字和下划线。
3. 双端 reads 的配对命名、singleton 处理和压缩格式符合所选 assembler 的要求。
4. probes 的 header 可由默认正则提取 locus；否则提前确定 `--regex`。
5. 用 `get-fastq-lengths` 抽查 reads 数量，用小型样本检查外部程序和日志路径。

### 6.1 组装 reads

SPAdes：

```bash
phyluce assembly assemblo-spades \
  --config samples.conf \
  --output assembly-spades \
  --subfolder split-adapter-quality-trimmed \
  --cores 12 \
  --memory 64
```

Velvet：

```bash
phyluce assembly assemblo-velvet \
  --config samples.conf \
  --output assembly-velvet \
  --kmer 35 \
  --subfolder split-adapter-quality-trimmed \
  --clean
```

ABySS：

```bash
phyluce assembly assemblo-abyss \
  --config samples.conf \
  --output assembly-abyss \
  --kmer 35 \
  --cores 12 \
  --subfolder split-adapter-quality-trimmed \
  --clean
```

样本配置示例：

```ini
[samples]
anas_platyrhynchos:/path/to/uce-clean/anas_platyrhynchos
gallus_gallus:/path/to/uce-clean/gallus_gallus
```

组装完成后，程序会尝试在输出目录下创建统一的 contigs 链接目录：

```text
assembly-output/
  contigs/
    anas_platyrhynchos.contigs.fasta
    gallus_gallus.contigs.fasta
```

注意：这些命令仍依赖 SPAdes、Velvet 或 ABySS 外部程序。

### 6.2 统计 FASTA/FASTQ 长度

FASTA：

```bash
phyluce assembly get-fasta-lengths --input sample.contigs.fasta
phyluce assembly get-fasta-lengths --input sample.contigs.fasta --csv
```

FASTQ：

```bash
phyluce assembly get-fastq-lengths --input raw-reads/
phyluce assembly get-fastq-lengths --input raw-reads/ --csv
```

### 6.3 contigs 匹配 probes

将 assembled contigs 与 UCE probes 做 LASTZ 匹配，并生成 `probe.matches.sqlite`：

```bash
phyluce assembly match-contigs-to-probes \
  --contigs assembly-spades/contigs \
  --probes uce-5k-probes.fasta \
  --output uce-search \
  --min-coverage 80 \
  --min-identity 80 \
  --csv uce-search/probe_match_results.csv
```

常用参数：

| 参数 | 说明 |
| --- | --- |
| `--contigs` | 包含 `*.fasta`, `*.fa`, `*.fna` contig 文件的目录 |
| `--probes` | UCE probe FASTA |
| `--output` | 输出目录，包含 `.lastz` 文件和 `probe.matches.sqlite` |
| `--min-coverage` | LASTZ 最低 coverage，默认 80 |
| `--min-identity` | LASTZ 最低 identity，默认 80 |
| `--regex` | 从 probe header 提取 locus 名称，默认 `^(uce-\d+)(?:_p\d+.*)` |
| `--dupefile` | probe 自比对 LASTZ 结果，用于移除疑似重复 probes |
| `--keep-duplicates` | 写出重复 locus/contig 信息 |
| `--csv` | 写出每个 taxon 的匹配摘要 |
| `--skip-alignment` | Rust 增加项：复用已有 `.lastz`，不调用 LASTZ |
| `--force` | Rust 增加项：输出目录已存在时删除并重建 |

输出：

```text
uce-search/
  sample1.contigs.lastz
  sample2.contigs.lastz
  probe.matches.sqlite
  probe_match_results.csv
```

匹配摘要不仅用于统计 UCE 回收量，也用于发现潜在旁系同源和重复命中。一个
locus 命中同一样本的多个 contig，或一个 contig 对应多个 locus 时，不应直接
当作多个独立 UCE 使用。应检查 `--dupefile`、重复命中输出和样本间回收量差异；
异常偏高可能来自重复序列或阈值过宽，异常偏低则可能反映数据质量、组装质量、
探针距离或命名问题。

### 6.4 生成 complete 或 incomplete matrix 的 locus 列表

配置文件示例：

```ini
[all]
alligator_mississippiensis
gallus_gallus
peromyscus_maniculatus
rana_sphenocephafa

[Excludes]
bad_sample
```

complete matrix：

```bash
phyluce assembly get-match-counts \
  --locus-db uce-search/probe.matches.sqlite \
  --taxon-list-config taxon-set.conf \
  --taxon-group all \
  --output taxon-set.complete.conf
```

incomplete matrix：

```bash
phyluce assembly get-match-counts \
  --locus-db uce-search/probe.matches.sqlite \
  --taxon-list-config taxon-set.conf \
  --taxon-group all \
  --incomplete-matrix \
  --output taxon-set.incomplete.conf
```

输出配置包含：

```ini
[Organisms]
...

[Loci]
...
```

这里的 complete matrix 要求所选 loci 在配置中的全部 taxa 都有匹配；
incomplete matrix 允许部分 taxa 缺失。先生成 incomplete 集合通常更便于评估
不同占比阈值会保留多少 loci，再根据研究设计选择最终矩阵。

### 6.5 从 contigs 提取 UCE FASTA

complete matrix：

```bash
phyluce assembly get-fastas-from-match-counts \
  --contigs assembly-spades/contigs \
  --locus-db uce-search/probe.matches.sqlite \
  --match-count-output taxon-set.complete.conf \
  --output taxon-set.complete.fasta
```

incomplete matrix：

```bash
phyluce assembly get-fastas-from-match-counts \
  --contigs assembly-spades/contigs \
  --locus-db uce-search/probe.matches.sqlite \
  --match-count-output taxon-set.incomplete.conf \
  --incomplete-matrix taxon-set.incomplete.missing \
  --output taxon-set.incomplete.fasta
```

该命令会：

1. 根据 `match_map` 找到每个 taxon 的 contig。
2. 按 `(+/-)` 方向必要时 reverse complement。
3. 将短 `N` 片段和 assembly 末端低覆盖小写碱基按原版逻辑处理。
4. 输出 monolithic FASTA，供下一步按 locus 比对。

### 6.6 拆分 monolithic FASTA

按 locus 拆分：

```bash
phyluce assembly explode-get-fastas-file \
  --input taxon-set.complete.fasta \
  --output unaligned-by-locus
```

按 taxon 拆分：

```bash
phyluce assembly explode-get-fastas-file \
  --input taxon-set.complete.fasta \
  --output unaligned-by-taxon \
  --by-taxon
```

### 6.7 阶段性质量检查

不要只以“命令成功退出”判断流程有效。每一步至少检查以下结果：

| 阶段 | 建议检查 | 常见异常信号 |
| --- | --- | --- |
| reads 与组装 | reads 数量、contig 数量和长度分布 | 某样本显著少于同批样本；大量极短 contig |
| probe 匹配 | 每个样本回收 loci 数、重复 locus/contig 数 | 单个样本回收量异常；重复命中过多 |
| match-count | taxa 名单、保留 loci 数、缺失分布 | 配置中的 taxon 未进入结果；阈值变化不合理 |
| FASTA 提取 | 序列条数、方向、header 和 missing 报告 | 空序列、异常短序列、header 无法拆分 |
| 比对与修剪 | 每个 locus 的 taxa 数、长度和模糊字符 | 修剪后大量 locus 消失或只剩很短片段 |
| 最终矩阵 | taxa 数、loci 数、总字符数和 charset 范围 | 拼接长度与 charset 末端不一致；样本全为缺失 |

建议把匹配摘要、alignment summary、taxon-locus counts 和最终 charset 与命令
日志一起归档。它们是定位数据损失发生在哪一阶段的最小证据链。

## 7. Alignment 工作流

MAFFT 是当前默认且推荐的比对程序。`seqcap-align` 可在比对后直接执行 edge
trimming；如需比较未修剪与已修剪结果，可对同一输入分别运行一次 `--no-trim`
和一次带修剪参数的命令。Gblocks/trimAl 可进一步处理内部不可靠区域，但更严格
的设置也可能删除具有系统发育信息的位点，因此阈值应结合数据和下游模型选择。

### 7.1 按 locus 比对 monolithic FASTA

```bash
phyluce align seqcap-align \
  --input taxon-set.complete.fasta \
  --output mafft-alignments \
  --taxa 4 \
  --no-trim
```

允许 incomplete matrix：

```bash
phyluce align seqcap-align \
  --input taxon-set.incomplete.fasta \
  --output mafft-alignments \
  --taxa 4 \
  --incomplete-matrix
```

启用 phyluce 原生 edge trimming：

```bash
phyluce align seqcap-align \
  --input taxon-set.complete.fasta \
  --output mafft-trimmed \
  --taxa 4 \
  --window 20 \
  --proportion 0.65 \
  --threshold 0.65 \
  --max-divergence 0.20 \
  --min-length 100 \
  --cores 8
```

注意：当前实现使用 MAFFT；`--aligner muscle` 尚未实现为 CLI 参数。
`seqcap-align` 无论是否启用修剪都写出 NEXUS。下一节的命令专用于已有的 FASTA
alignment 目录，不能直接读取这里生成的 NEXUS；如需这样串联，应先转换为 FASTA。

### 7.2 修剪已有 FASTA alignments

phyluce 原生三阶段 edge trimming：

```bash
phyluce align get-trimmed-alignments-from-untrimmed \
  --alignments mafft-alignments \
  --output mafft-edge-trim \
  --window 20 \
  --proportion 0.65 \
  --threshold 0.65 \
  --max-divergence 0.20 \
  --min-length 100 \
  --cores 8
```

Gblocks：

```bash
phyluce align get-gblocks-trimmed-alignments-from-untrimmed \
  --alignments mafft-alignments \
  --output mafft-gblocks \
  --input-format fasta \
  --output-format nexus \
  --cores 8
```

trimAl：

```bash
phyluce align get-trimal-trimmed-alignments-from-untrimmed \
  --alignments mafft-alignments \
  --output mafft-trimal \
  --input-format fasta \
  --output-format nexus \
  --cores 8
```

这些命令按 `--cores` 并行处理独立 locus；Gblocks 和 trimAl 仍依赖外部程序。

修剪后运行 `get-align-summary-data`，对比输入与输出的 locus 数、平均长度和
taxa 数。如果大量 locus 被删除或长度骤降，应先复核参数，不要直接进入拼接。

### 7.3 添加缺失数据占位符

```bash
phyluce align add-missing-data-designators \
  --alignments mafft-alignments \
  --output mafft-with-missing \
  --match-count-output taxon-set.incomplete.conf \
  --incomplete-matrix taxon-set.incomplete.missing \
  --input-format fasta
```

默认缺失字符为 `?`：

```bash
--missing-character ?
```

### 7.4 移除空 taxon

```bash
phyluce align remove-empty-taxa \
  --alignments mafft-with-missing \
  --output mafft-no-empty \
  --input-format nexus \
  --output-format nexus \
  --cores 8
```

### 7.5 alignment 格式转换

```bash
phyluce align convert-one-align-to-another \
  --alignments mafft-alignments \
  --output mafft-nexus \
  --input-format fasta \
  --output-format nexus \
  --cores 8
```

输入支持 FASTA、NEXUS、PHYLIP（含 relaxed/sequential）、CLUSTAL、EMBOSS
和 Stockholm；该转换命令当前仅输出 FASTA 或 NEXUS，不支持的格式会明确报错。

### 7.6 转换 IUPAC degenerate bases

```bash
phyluce align convert-degen-bases \
  --alignments mafft-degen-bases \
  --output mafft-degen-bases-converted \
  --input-format fasta \
  --output-format nexus \
  --cores 8
```

该命令将 degenerate IUPAC 碱基转换为 `N`。

### 7.7 筛选和提取 taxa

从 alignments 中排除 taxon：

```bash
phyluce align extract-taxa-from-alignments \
  --alignments mafft-gblocks-clean \
  --output no-gallus \
  --input-format nexus \
  --output-format nexus \
  --exclude gallus_gallus
```

只保留指定 taxa：

```bash
phyluce align extract-taxa-from-alignments \
  --alignments mafft-gblocks-clean \
  --output keep-two \
  --input-format nexus \
  --output-format nexus \
  --include gallus_gallus peromyscus_maniculatus
```

提取某个 taxon 的 FASTA：

```bash
phyluce align extract-taxon-fasta-from-alignments \
  --alignments mafft-gblocks-clean \
  --taxon gallus_gallus \
  --output gallus.fasta \
  --input-format nexus
```

### 7.8 根据条件过滤 alignment

```bash
phyluce align filter-alignments \
  --alignments mafft-gblocks-clean \
  --output filtered-alignments \
  --input-format nexus \
  --containing-data-for gallus_gallus \
  --min-length 600 \
  --min-taxa 3 \
  --cores 8
```

按最小 taxon 占比保留 loci：

```bash
phyluce align get-only-loci-with-min-taxa \
  --alignments mafft-gblocks-clean \
  --taxa 4 \
  --percent 0.75 \
  --output mafft-gblocks-clean-75p \
  --input-format nexus \
  --cores 8
```

`--percent 0.75` 表示每个保留 locus 至少包含 `--taxa` 所定义总 taxon 数的
75%，不是“每个 taxon 至少出现在 75% 的 loci 中”。例如 `--taxa 20
--percent 0.75` 要求每个 locus 至少有 15 个 taxa；它不能保证任一特定 taxon
在最终矩阵中的覆盖率。后者应通过 `get-taxon-locus-counts-in-alignments` 单独检查。

### 7.9 清理序列名并拼接 alignments

从 UCE monolithic FASTA 生成的序列名可能包含 locus 后缀。下游分析通常只需要
taxon 名，拼接前可先清理：

```bash
phyluce align remove-locus-name-from-files \
  --alignments mafft-gblocks-clean \
  --output cleaned-names \
  --input-format nexus \
  --output-format nexus \
  --cores 8
```

先抽查输出 header，确认清理后没有产生重复 taxon 名，再进行拼接。

输出 NEXUS：

```bash
phyluce align concatenate-alignments \
  --alignments cleaned-names \
  --output concat-nexus \
  --nexus
```

输出 PHYLIP 和 charset：

```bash
phyluce align concatenate-alignments \
  --alignments cleaned-names \
  --output concat-phylip \
  --phylip
```

拼接过程分两遍读取 locus，并使用自动清理的磁盘暂存矩阵；峰值内存主要由单个
locus 决定，而不是全部输入和完整拼接矩阵之和。

PHYLIP 输出同时生成 charset 信息，记录每个 locus 在拼接矩阵中的字符范围，可
作为 IQ-TREE、RAxML 等分区配置的依据；是否需要转换格式取决于下游软件版本和
参数。运行下游程序前，应确认所有 charset 区间连续、不重叠，最后一个区间终点
等于拼接序列总长度。

### 7.10 拆分 concatenated NEXUS

```bash
phyluce align split-concat-nexus-to-loci \
  --nexus concat-nexus/concat-nexus.nexus \
  --output split-loci \
  --output-format nexus
```

### 7.11 统计 alignment

alignment summary：

```bash
phyluce align get-align-summary-data \
  --alignments mafft-gblocks-clean \
  --input-format nexus \
  --cores 8 \
  --output-stats align-summary.csv
```

alignment 文件按 `--cores` 并行解析和统计，最终输出仍按文件名排序。

informative sites：

```bash
phyluce align get-informative-sites \
  --alignments mafft-gblocks-clean \
  --input-format nexus \
  --output informative-sites.csv
```

taxon-locus counts：

```bash
phyluce align get-taxon-locus-counts-in-alignments \
  --alignments mafft-gblocks-clean \
  --input-format nexus \
  --output taxon-locus-counts.csv
```

incomplete matrix 估计：

```bash
phyluce align get-incomplete-matrix-estimates \
  --db uce-search/probe.matches.sqlite \
  --min 0.5 \
  --max 1.0 \
  --step 0.05
```

### 7.12 RY recoding

RY：

```bash
phyluce align get-ry-recoded-alignments \
  --alignments mafft-gblocks-clean \
  --output mafft-ry \
  --input-format nexus
```

二进制编码：

```bash
phyluce align get-ry-recoded-alignments \
  --alignments mafft-gblocks-clean \
  --output mafft-ry-binary \
  --input-format nexus \
  --binary
```

### 7.13 其他 alignment 命令

按配置移动 alignments：

```bash
phyluce align move-align-by-conf-file \
  --conf selected-loci.conf \
  --alignments mafft-gblocks-clean \
  --output selected-alignments \
  --sections all \
  --extension nex
```

随机抽样并拼接：

```bash
phyluce align randomly-sample-and-concatenate \
  --alignments mafft-gblocks-clean \
  --output sampled-concat \
  --sample-size 100 \
  --replicates 10
```

使用 RAxML reduction：

```bash
phyluce align reduce-alignments-with-raxml \
  --alignments phylip-alignments \
  --output reduced-alignments \
  --input-format phylip-relaxed
```

移除序列名中的 locus 名称：

```bash
phyluce align remove-locus-name-from-files \
  --alignments mafft-gblocks-clean \
  --output cleaned-names \
  --input-format nexus \
  --output-format nexus
```

筛查异常碱基：

```bash
phyluce align screen-alignments-for-problems \
  --alignments mafft-gblocks-clean \
  --output screened-alignments \
  --input-format nexus
```

生成 smilogram 数据：

```bash
phyluce align get-smilogram-from-alignments \
  --alignments mafft-gblocks-clean \
  --output-file smilogram.csv \
  --output-missing smilogram-missing.csv \
  --output-database smilogram.sqlite \
  --input-format nexus
```

为 PAML 格式化 concatenated PHYLIP：

```bash
phyluce align format-concatenated-phylip-for-paml \
  --phylip-alignment concat.phylip \
  --config charsets.conf \
  --output paml-ready.phylip
```

## 8. Probe 命令

### 8.1 从 LASTZ 结果生成 BED

probe-level BED：

```bash
phyluce probe get-probe-bed-from-lastz-files \
  --alignments lastz-output \
  --output probe-bed-output
```

locus-level BED：

```bash
phyluce probe get-locus-bed-from-lastz-files \
  --alignments lastz-output \
  --output locus-bed-output \
  --regex '^(uce-\d+)(?:_p\d+.*)'
```

### 8.2 根据配置过滤 probes

```bash
phyluce probe remove-overlapping-probes-given-config \
  --probes probes.fasta \
  --config keep.conf \
  --output filtered-probes.fasta
```

### 8.3 筛选 tiled probes 子集

```bash
phyluce probe get-subsets-of-tiled-probes \
  --probes tiled-probes.fasta \
  --taxa alligator gallus \
  --output subset-probes.fasta
```

### 8.4 多 FASTA / multi-merge SQLite 表

生成 multi-fasta table：

```bash
phyluce probe get-multi-fasta-table \
  --fastas fastas/ \
  --output multi-fasta.sqlite \
  --base-taxon alligator
```

查询 multi-fasta table：

```bash
phyluce probe query-multi-fasta-table \
  --db multi-fasta.sqlite \
  --base-taxon alligator \
  --specific-counts 3
```

生成 multi-merge table：

```bash
phyluce probe get-multi-merge-table \
  --conf genomes.conf \
  --output multi-merge.sqlite \
  --base-taxon alligator
```

查询 multi-merge table：

```bash
phyluce probe query-multi-merge-table \
  --db multi-merge.sqlite \
  --base-taxon alligator \
  --specific-counts 3
```

### 8.5 根据距离筛选 loci

```bash
phyluce probe get-screened-loci-by-proximity \
  --input probes.fasta \
  --output screened-loci.fasta \
  --distance 10000
```

注意：原版用随机方式解决同一 cluster 内的 tie；Rust 版本选择最小 locus id，结果可重复。

### 8.6 移除重复 hits

```bash
phyluce probe remove-duplicate-hits-from-probes-using-lastz \
  --fasta probes.fasta \
  --lastz probes-self.lastz \
  --probe-prefix uce- \
  --long
```

可选输出 BED：

```bash
--probe-bed probes.bed --locus-bed loci.bed
```

### 8.7 设计 tiled probes

单输入：

```bash
phyluce probe get-tiled-probes \
  --input loci.fasta \
  --output probes.fasta \
  --probe-prefix uce- \
  --designer faircloth \
  --design test-design \
  --probe-length 120 \
  --tiling-density 2 \
  --overlap middle \
  --probe-bed probes.bed \
  --locus-bed loci.bed
```

多输入：

```bash
phyluce probe get-tiled-probe-from-multiple-inputs \
  --fastas loci-by-taxon \
  --multi-fasta-output multi-fasta.conf \
  --output tiled-probes.fasta \
  --probe-prefix uce- \
  --designer faircloth \
  --design test-design \
  --probe-length 120 \
  --tiling-density 2
```

注意：`--two-probes` 的 tie 处理在 Rust 中是确定性的，不再使用随机选择。

### 8.8 从 probes 重建 UCE

```bash
phyluce probe reconstruct-uce-from-probe \
  --input probes.fasta \
  --output reconstructed-uces.fasta
```

如果一个 locus 有多个 probe，默认使用配置中的 MAFFT：

```bash
phyluce probe reconstruct-uce-from-probe \
  --input probes.fasta \
  --output reconstructed-uces.fasta
```

可通过 `--mafft-binary /path/to/mafft` 覆盖默认路径。需要复现原版
MUSCLE 3/Clustal alignment 路径时，显式传入
`--muscle-binary /path/to/muscle`。

### 8.9 2bit / BED / genome sequence 工具

从 BED 提取 genome sequence：

```bash
phyluce probe get-genome-sequences-from-bed \
  --bed loci.bed \
  --twobit genome.2bit \
  --output loci.fasta \
  --filter-mask 0.25 \
  --max-n 0
```

过滤 masked loci：

```bash
phyluce probe strip-masked-loci-from-set \
  --bed loci.bed \
  --twobit genome.2bit \
  --output unmasked-loci.fasta \
  --filter-mask 0.25 \
  --max-n 0 \
  --min-length 100
```

根据 LASTZ 结果从多个 genome 中切片：

```bash
phyluce probe slice-sequence-from-genomes \
  --conf genomes.conf \
  --lastz multi.lastz \
  --output sliced-fastas \
  --flank 500
```

配置示例：

```ini
[chromos]
alligator:/path/to/alligator.2bit
gallus:/path/to/gallus.2bit

[scaffolds]
taxon_x:/path/to/taxon_x.2bit
```

必须且只能选择一个：

```bash
--flank 500
--probes 3
```

### 8.10 LASTZ probe 命令

简单 LASTZ：

```bash
phyluce probe easy-lastz \
  --target target.fasta \
  --query probes.fasta \
  --output result.lastz \
  --identity 92.5 \
  --coverage 83
```

多 genome LASTZ 并写 SQLite：

```bash
phyluce probe run-multiple-lastzs-sqlite \
  --db multi-lastz.sqlite \
  --output lastz-output \
  --probefile probes.fasta \
  --genome-base-path /path/to/genomes \
  --chromolist alligator gallus \
  --cores 8 \
  --coverage 83 \
  --identity 92.5
```

`--cores` 控制全局并发 LASTZ 进程数。所有 genome 共用一个有界任务队列：
染色体 genome 按 `.2bit` 内的序列拆分，scaffold genome 边解码边生成约 10 Mbp
临时 FASTA 分块。各分块可以乱序完成，但会按目标顺序合并；完成的 genome 会
立即流式生成 `.clean` 文件并由主线程逐行写入 SQLite 事务。

## 9. Utilities 命令

从 FASTA header 生成 BED：

```bash
phyluce utilities get-bed-from-fasta \
  --input loci.fasta \
  --output loci.bed \
  --locus-prefix uce-
```

根据 FASTA 过滤 BED：

```bash
phyluce utilities filter-bed-by-fasta \
  --bed loci.bed \
  --fasta keep.fasta \
  --output filtered.bed
```

批量替换 links：

```bash
phyluce utilities replace-many-links \
  --indir links-in \
  --oldpath /old/base \
  --newpath /new/base \
  --outdir links-out
```

合并 reads：

```bash
phyluce utilities combine-reads \
  --config reads.conf \
  --output combined-reads \
  --subfolder split-adapter-quality-trimmed
```

合并多个 gzip：

```bash
phyluce utilities merge-multiple-gzip-files \
  --config samples.conf \
  --output merged \
  --section samples
```

合并 NextSeq gzip：

```bash
phyluce utilities merge-next-seq-gzip-files \
  --input nextseq-run \
  --config samples.conf \
  --output merged-nextseq \
  --section samples
```

拆分 mixed FASTA reads：

```bash
phyluce utilities unmix-fasta-reads \
  --mixed-reads mixed.fasta \
  --out-r1 reads_R1.fasta \
  --out-r2 reads_R2.fasta \
  --out-r-singleton reads_singleton.fasta
```

用 seqtk 抽样 reads：

```bash
phyluce utilities sample-reads-from-files \
  --conf samples.conf \
  --output sampled-reads
```

注意：`sample-reads-from-files` 仍依赖外部 `seqtk`。

## 10. NCBI 命令

将 FASTA 拆成 NCBI 可接受的块：

```bash
phyluce ncbi chunk-fasta-for-ncbi \
  --input all-uces.fasta \
  --chunk-size 10000 \
  --output-prefix split \
  --output-suffix fsa
```

准备 UCE alignment FASTA 供 NCBI 提交：

```bash
phyluce ncbi prep-uce-align-files-for-ncbi \
  --alignments alignments \
  --conf ncbi-prep.conf \
  --output ncbi-output \
  --input-format nexus
```

注意：原版该命令在现代 Biopython 下可能因 `Bio.Alphabet` 删除而无法导入；Rust 版本按预期行为实现。

## 11. Genetrees 命令

重命名 tree leaves：

```bash
phyluce genetrees rename-tree-leaves \
  --input tree.tre \
  --config names.conf \
  --section standard \
  --output renamed.tre
```

tree topology 计数：

```bash
phyluce genetrees get-tree-counts \
  --trees gene-trees \
  --locus-support-output locus-support.csv \
  --root outgroup_taxon \
  --extension tre
```

平均 bootstrap replicate support：

```bash
phyluce genetrees get-mean-bootrep-support \
  --trees bootstrap-trees \
  --config taxon-map.conf
```

输出文件名固定为当前目录下的 `outfile.csv`，与原版一致。

生成 multilocus bootstrap count：

```bash
phyluce genetrees generate-multilocus-bootstrap-count \
  --alignments alignments \
  --bootstrap-replicates replicates.txt \
  --bootstrap-counts counts.txt \
  --bootreps 100
```

排序 multilocus bootstraps：

```bash
phyluce genetrees sort-multilocus-bootstraps \
  --input bootstrap-input \
  --bootstrap-replicates replicates.txt \
  --output sorted-bootstrap
```

注意：Rust 版本的 bootstrap replicate 文件使用纯文本格式，不使用 Python `pickle`。

## 12. Workflow 命令

Rust 版本保留原版 `phyluce_workflow` 的入口，用于调用 Snakemake workflow：

```bash
phyluce workflow \
  --config mapping.config.yaml \
  --output workflow-output \
  --workflow mapping \
  --cores 8
```

支持 workflow：

```text
mapping
correction
phasing
```

dry run：

```bash
phyluce workflow \
  --config mapping.config.yaml \
  --output workflow-output \
  --workflow mapping \
  --cores 8 \
  --dryrun
```

注意：该命令仍需要 Snakemake 和对应外部程序。

## 13. I/O 和诊断命令

检查 FASTA：

```bash
phyluce io validate-fasta --input probes.fasta
```

如果 FASTA 格式正确，会输出：

```text
OK: probes.fasta is well-formed FASTA
```

检查外部程序：

```bash
phyluce external check --program binaries --binary mafft
```

输出包括解析路径、退出码和 `--version` 输出。

## 14. 已知差异

| 命令/功能 | 差异 |
| --- | --- |
| `assembly match-contigs-to-barcodes` | 不执行 BOLD 网络查询；本地 LASTZ slicing 请使用 `--no-bold` |
| `assembly match-contigs-to-probes` | 新增 `--skip-alignment` 和 `--force`，用于 CI 和非交互运行 |
| `align seqcap-align` | 当前使用 MAFFT；原版支持 MAFFT/MUSCLE 选择 |
| alignment 输入 | 支持 FASTA、NEXUS、PHYLIP（含 relaxed/sequential）、CLUSTAL、EMBOSS 和 Stockholm；各命令的输出格式仍以帮助信息和明确报错为准 |
| `align randomly-sample-and-concatenate` | 使用 seeded PRNG 思路，避免原版随机行为不可复现 |
| `align get-smilogram-from-alignments` | major-allele tie 使用确定性规则 |
| `probe get-screened-loci-by-proximity` | cluster tie 保留最小 locus id，而非随机选择 |
| `probe get-tiled-probes` / `get-tiled-probe-from-multiple-inputs` | `--two-probes` tie 处理为确定性 |
| `probe reconstruct-uce-from-probe` | 默认使用 MAFFT；可通过 `--muscle-binary` 显式使用原版 MUSCLE 3/Clustal 路径 |
| `genetrees generate-multilocus-bootstrap-count` | 使用纯文本 replicate 格式，不使用 Python pickle |
| `assembly get-match-counts` | 尚未移植原版 `--optimize` 随机优化路径 |
| `genetrees rename-tree-leaves` | 尚未实现 `--reroot`；部分 genetree 命令仅接受 Newick 输入 |
| `ncbi prep-uce-align-files-for-ncbi` | Rust 版按预期行为实现，不复现现代 Biopython 下原版导入失败 |

## 15. 故障排查

### 15.1 找不到外部程序

先查看配置：

```bash
phyluce config inspect
phyluce config which --program binaries --binary lastz
```

再检查程序：

```bash
phyluce external check --program binaries --binary lastz
```

如果路径仍指向错误 conda 环境，请更新 `~/.phyluce.conf`。

### 15.2 输出目录已存在

部分原版脚本会交互式询问是否删除输出目录。Rust 版本通常不做交互删除。对于支持的命令，可以使用：

```bash
--force
```

或先手动换一个输出目录。

### 15.3 LASTZ 不可用但已有 `.lastz` 文件

`assembly match-contigs-to-probes` 支持复用已有 LASTZ 输出：

```bash
phyluce assembly match-contigs-to-probes \
  --contigs contigs \
  --probes probes.fasta \
  --output uce-search \
  --skip-alignment
```

### 15.4 FASTA header 无法解析

先确认 FASTA 格式：

```bash
phyluce io validate-fasta --input input.fasta
```

对于 probes，默认 locus 提取正则为：

```text
^(uce-\d+)(?:_p\d+.*)
```

如果 probe 命名不同，请显式传入 `--regex` 或对应命令中的 `--probe-regex`。

### 15.5 与原版输出不完全一致

优先检查：

1. 外部工具版本是否一致，尤其是 LASTZ、MAFFT、Gblocks、trimAl。
2. 输入文件排序是否一致。
3. 是否使用了 Rust 版本的确定性 tie 处理。
4. 是否启用了 `--skip-alignment`、`--force` 等 Rust 增加项。

## 16. 兼容性测试

独立仓库可直接运行内置 fixture：

```bash
python3 tests/compat/run_fixtures.py target/debug/phyluce
```

与原版 Python 做完整动态对照时，指定原版源码目录：

```bash
PHYLUCE_PYTHON_REPO=/path/to/phyluce \
  python3 tests/compat/run_all.py target/debug/phyluce
```

完整动态对照会调用原版依赖以及 MAFFT、LASTZ、RAxML 等外部程序；缺少这些
依赖时应使用 `run_fixtures.py`。

## 17. 端到端最小示例

以下示例从已组装 contigs 开始，构建一个至少包含 75% taxa 的矩阵。假设
`taxon-sets/taxon-set.conf` 的 `[all]` 中有 4 个 taxa；真实分析必须把 `--taxa 4`
改为配置中的总数。每一步使用独立目录，并把日志保存在 `logs/`。

先检查输入和外部程序：

```bash
mkdir -p logs stats taxon-sets alignments matrices
phyluce --version
phyluce config inspect
phyluce external check --program binaries --binary lastz
phyluce external check --program binaries --binary mafft
phyluce io validate-fasta --input probes/uce-5k-probes.fasta
```

运行主流程：

```bash
# 1. contigs 与 probes 匹配
phyluce --log-path logs/01-match assembly match-contigs-to-probes \
  --contigs assembly/contigs \
  --probes probes/uce-5k-probes.fasta \
  --output uce-search \
  --csv uce-search/probe_match_results.csv

# 2. 生成允许缺失的 locus 集合
phyluce --log-path logs/02-counts assembly get-match-counts \
  --locus-db uce-search/probe.matches.sqlite \
  --taxon-list-config taxon-sets/taxon-set.conf \
  --taxon-group all \
  --incomplete-matrix \
  --output taxon-sets/all.incomplete.conf

# 3. 提取并重命名 UCE 序列
phyluce --log-path logs/03-extract assembly get-fastas-from-match-counts \
  --contigs assembly/contigs \
  --locus-db uce-search/probe.matches.sqlite \
  --match-count-output taxon-sets/all.incomplete.conf \
  --incomplete-matrix taxon-sets/all.missing.conf \
  --output taxon-sets/all.incomplete.fasta

# 4. MAFFT 比对并直接执行 edge trimming；输出为 NEXUS
phyluce --log-path logs/04-align align seqcap-align \
  --input taxon-sets/all.incomplete.fasta \
  --output alignments/trimmed \
  --taxa 4 \
  --incomplete-matrix \
  --cores 8

# 5. 汇总修剪结果
phyluce --log-path logs/05-summary align get-align-summary-data \
  --alignments alignments/trimmed \
  --input-format nexus \
  --cores 8 \
  --output-stats stats/trimmed-alignments.csv

# 6. 保留每个 locus 至少含 75% taxa 的 alignment
phyluce --log-path logs/06-filter align get-only-loci-with-min-taxa \
  --alignments alignments/trimmed \
  --taxa 4 \
  --percent 0.75 \
  --output alignments/filtered-75p \
  --input-format nexus \
  --cores 8

# 7. 检查每个 taxon 在多少 loci 中有数据
phyluce align get-taxon-locus-counts-in-alignments \
  --alignments alignments/filtered-75p \
  --input-format nexus \
  --output stats/taxon-locus-counts-75p.csv

# 8. 移除序列名中的 locus 后缀
phyluce --log-path logs/08-clean-names align remove-locus-name-from-files \
  --alignments alignments/filtered-75p \
  --output alignments/final-75p \
  --input-format nexus \
  --output-format nexus \
  --cores 8

# 9. 输出 PHYLIP 与 charset/partition 文件
phyluce --log-path logs/09-concat align concatenate-alignments \
  --alignments alignments/final-75p \
  --output matrices/all-75p-phylip \
  --phylip
```

完成后至少核对：`probe_match_results.csv` 中各样本回收量没有异常离群值；
`trimmed-alignments.csv` 中没有大批极短 loci；每个最终 alignment 满足 75%
阈值；charset 的最后坐标等于拼接长度；关键 taxon 在
`taxon-locus-counts-75p.csv` 中没有大面积缺失。需要 NEXUS 时，对同一个
`alignments/final-75p` 再运行一次 `concatenate-alignments --nexus`，并使用不同
输出目录。

## 18. 版本状态

当前 workspace 版本为 `0.1.0`，Rust edition 为 2021。建议在正式数据分析中记录：

```bash
phyluce --version
phyluce config inspect
```

并保存每一步命令的 `--log-path` 日志，以便复现。

## 19. 可复现性与资源设置

正式分析应保存以下信息：

```bash
phyluce --version
phyluce config inspect
git rev-parse HEAD               # 从源码构建时记录
mafft --version
lastz --version
```

同时归档：

- 样本与 taxon 配置、probe FASTA 的版本或校验值。
- 每一步完整命令、日志、输入和输出目录名称。
- MAFFT、LASTZ、assembler、trimming 和系统发育软件的版本。
- filtering、occupancy、identity、coverage 和 trimming 参数。
- 最终 taxa/loci 数、alignment 统计、charset/partition 和缺失数据摘要。

`--cores` 控制支持并行的命令，但最优值不一定等于机器全部逻辑核心。MAFFT、
LASTZ、组装器和文件级并行可能同时消耗 CPU、内存与文件句柄；先用少量样本观察
资源占用，再增加核心数。多 genome LASTZ 任务使用全局有界任务队列，拼接使用
两遍读取和磁盘暂存，因此还应确保输入文件在运行期间不被修改，并为系统临时
目录（通常由 `TMPDIR` 决定）保留足够空间。

## 20. 引用、许可证与问题报告

使用本项目进行研究时，应记录 Rust CLI 版本或提交号，并根据实际工作流引用
原版 PHYLUCE 方法及所调用的外部软件。原版 PHYLUCE 建议引用：

> Faircloth, B. C. 2016. PHYLUCE is a software package for the analysis of
> conserved genomic loci. *Bioinformatics* 32:786-788.
> https://doi.org/10.1093/bioinformatics/btv646

探针设计、目标富集方法和特定生物类群可能还需要其他论文；以
[原版引用说明](https://phyluce.readthedocs.io/en/latest/citing.html)为准。Rust
重构本身当前没有单独的论文或 DOI，不应虚构引用信息。

本项目采用 BSD 3-Clause License，完整条款见 [LICENSE](../LICENSE)。分发源码或
二进制时应遵守其中的版权声明、条件和免责声明。

报告问题时，请在当前仓库的
[GitHub Issues](https://github.com/GUIBA-EX/phyluce_rust/issues)提供：

1. `phyluce --version`、操作系统和 CPU 架构。
2. 完整命令、最小可复现输入和使用的配置片段；删除敏感路径。
3. 对应日志、实际结果、预期结果和退出状态。
4. 涉及的外部程序及其版本。
5. 问题是否可在单核心和小型 fixture 上复现。

## 21. 参考资料

- [PHYLUCE 官方文档](https://phyluce.readthedocs.io/en/latest/index.html)：原版项目
  的整体入口和命令索引。
- [Purpose](https://phyluce.readthedocs.io/en/latest/purpose.html)：UCE 工作流的
  目标和适用范围。
- [Installation](https://phyluce.readthedocs.io/en/latest/installation.html)：原版
  Python 环境及外部依赖背景；不要直接作为 Rust 安装步骤使用。
- [Tutorial I](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html)：
  从 reads、组装、UCE 提取到矩阵准备的经典流程。
- [Citing](https://phyluce.readthedocs.io/en/latest/citing.html)：PHYLUCE 与 UCE
  方法的引用建议。
- [License](https://phyluce.readthedocs.io/en/latest/license.html)：原版软件和文档
  的许可信息。本项目的实际许可仍以本仓库 `LICENSE` 为准。
