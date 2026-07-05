# phyluce Rust CLI 中文用户手册

本文档面向已经熟悉原版 `phyluce` 的用户，说明 Rust 版本如何安装、配置和使用。Rust 版本的目标是模仿原版脚本的工作方式，同时提供一个更统一的命令入口。

当前 Rust CLI 处于早期预览阶段。多数命令已经按原版脚本行为实现，但仍有少数命令依赖外部程序，或与原版存在有意的兼容差异。正式分析前，建议先用小数据或现有 fixture 运行一遍完整流程。

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

Rust 版本会解析默认配置和用户配置，并展开 `$CONDA`、`$WORKFLOWS`。

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

Rust 二进制支持部分旧脚本名兼容：如果把 `phyluce` 复制或软链接为旧脚本名，它会自动映射到新的分组命令。

例如：

```bash
ln -s target/release/phyluce phyluce_align_convert_degen_bases
./phyluce_align_convert_degen_bases --alignments in --output out
```

等价于：

```bash
phyluce align convert-degen-bases --alignments in --output out
```

当前已映射的旧脚本名：

```text
phyluce_align_convert_degen_bases
phyluce_align_explode_alignments
phyluce_align_extract_taxon_fasta_from_alignments
phyluce_align_format_concatenated_phylip_for_paml
phyluce_align_get_incomplete_matrix_estimates
phyluce_align_get_only_loci_with_min_taxa
phyluce_align_get_taxon_locus_counts_in_alignments
phyluce_align_move_align_by_conf_file
phyluce_align_randomly_sample_and_concatenate
phyluce_align_reduce_alignments_with_raxml
phyluce_align_remove_locus_name_from_files
phyluce_align_screen_alignments_for_problems
phyluce_align_get_smilogram_from_alignments
phyluce_assembly_screen_probes_for_dupes
phyluce_assembly_extract_contigs_to_barcodes
phyluce_assembly_match_contigs_to_barcodes
```

未列出的命令请使用新的分组式 CLI。

## 6. UCE 主流程

典型 UCE phylogenomics 流程如下：

```text
raw reads
  -> assembly
  -> contigs/*.contigs.fasta
  -> contigs 与 probes 做 LASTZ 匹配
  -> probe.matches.sqlite
  -> match-count config
  -> monolithic UCE FASTA
  -> 按 locus 比对
  -> trimming / filtering / concatenation / summary
```

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

## 7. Alignment 工作流

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
  --min-length 100
```

注意：当前实现使用 MAFFT；`--aligner muscle` 尚未实现为 CLI 参数。

### 7.2 修剪已有 alignments

phyluce 原生三阶段 edge trimming：

```bash
phyluce align get-trimmed-alignments-from-untrimmed \
  --alignments mafft-alignments \
  --output mafft-edge-trim \
  --window 20 \
  --proportion 0.65 \
  --threshold 0.65 \
  --max-divergence 0.20 \
  --min-length 100
```

Gblocks：

```bash
phyluce align get-gblocks-trimmed-alignments-from-untrimmed \
  --alignments mafft-alignments \
  --output mafft-gblocks \
  --input-format fasta \
  --output-format nexus
```

trimAl：

```bash
phyluce align get-trimal-trimmed-alignments-from-untrimmed \
  --alignments mafft-alignments \
  --output mafft-trimal \
  --input-format fasta \
  --output-format nexus
```

Gblocks 和 trimAl 仍依赖外部程序。

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
  --output-format nexus
```

### 7.5 alignment 格式转换

```bash
phyluce align convert-one-align-to-another \
  --alignments mafft-alignments \
  --output mafft-nexus \
  --input-format fasta \
  --output-format nexus
```

当前重点支持 FASTA/NEXUS 兼容路径。其他格式请先用小数据验证。

### 7.6 转换 IUPAC degenerate bases

```bash
phyluce align convert-degen-bases \
  --alignments mafft-degen-bases \
  --output mafft-degen-bases-converted \
  --input-format fasta \
  --output-format nexus
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
  --min-taxa 3
```

按最小 taxon 占比保留 loci：

```bash
phyluce align get-only-loci-with-min-taxa \
  --alignments mafft-gblocks-clean \
  --taxa 4 \
  --percent 0.75 \
  --output mafft-gblocks-clean-75p \
  --input-format nexus
```

### 7.9 拼接 alignments

输出 NEXUS：

```bash
phyluce align concatenate-alignments \
  --alignments mafft-gblocks-clean \
  --output concat-nexus \
  --nexus
```

输出 PHYLIP 和 charset：

```bash
phyluce align concatenate-alignments \
  --alignments mafft-gblocks-clean \
  --output concat-phylip \
  --phylip
```

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
  --output-stats align-summary.csv
```

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

如果一个 locus 有多个 probe，需要 MAFFT：

```bash
phyluce probe reconstruct-uce-from-probe \
  --input probes.fasta \
  --output reconstructed-uces.fasta \
  --mafft-binary /path/to/mafft
```

注意：原版多 probe locus 使用 MUSCLE/Clustal 路径；Rust 版本使用 MAFFT。

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
  --coverage 83 \
  --identity 92.5
```

注意：当前 Rust 版本接受 `--cores`，但该实现未按原版 multiprocessing/chunking 并行化。

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
| `align convert-one-align-to-another` | 当前主要覆盖 FASTA/NEXUS 兼容路径 |
| `align randomly-sample-and-concatenate` | 使用 seeded PRNG 思路，避免原版随机行为不可复现 |
| `align get-smilogram-from-alignments` | major-allele tie 使用确定性规则 |
| `probe get-screened-loci-by-proximity` | cluster tie 保留最小 locus id，而非随机选择 |
| `probe get-tiled-probes` / `get-tiled-probe-from-multiple-inputs` | `--two-probes` tie 处理为确定性 |
| `probe reconstruct-uce-from-probe` | 多 probe locus 使用 MAFFT，而非原版 MUSCLE/Clustal 路径 |
| `probe run-multiple-lastzs-sqlite` | `--cores` 已接受但未按原版 chunked multiprocessing 并行化 |
| `genetrees generate-multilocus-bootstrap-count` | 使用纯文本 replicate 格式，不使用 Python pickle |
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

仓库包含兼容性测试脚本：

```bash
python3 tests/compat/run_all.py target/debug/phyluce
```

也可以单独运行某个测试：

```bash
python3 tests/compat/compare_get_fasta_lengths.py target/debug/phyluce
python3 tests/compat/compare_match_contigs_to_probes.py target/debug/phyluce
python3 tests/compat/compare_seqcap_align.py target/debug/phyluce
```

注意：完整测试会调用 MAFFT、LASTZ、RAxML 等外部程序。若当前环境没有这些程序，请只运行不依赖外部程序的子集，或在完整 conda 环境中运行。

## 17. 推荐最小示例

以下命令展示一个最小 UCE 主流程：

```bash
# 1. contigs vs probes
phyluce assembly match-contigs-to-probes \
  --contigs assembly/contigs \
  --probes uce-5k-probes.fasta \
  --output uce-search \
  --csv uce-search/probe_match_results.csv

# 2. complete matrix loci
phyluce assembly get-match-counts \
  --locus-db uce-search/probe.matches.sqlite \
  --taxon-list-config taxon-set.conf \
  --taxon-group all \
  --output taxon-set.complete.conf

# 3. extract monolithic FASTA
phyluce assembly get-fastas-from-match-counts \
  --contigs assembly/contigs \
  --locus-db uce-search/probe.matches.sqlite \
  --match-count-output taxon-set.complete.conf \
  --output taxon-set.complete.fasta

# 4. align by locus
phyluce align seqcap-align \
  --input taxon-set.complete.fasta \
  --output mafft-no-trim \
  --taxa 4 \
  --no-trim

# 5. trim existing alignments
phyluce align get-trimmed-alignments-from-untrimmed \
  --alignments mafft-no-trim \
  --output mafft-edge-trim

# 6. concatenate
phyluce align concatenate-alignments \
  --alignments mafft-edge-trim \
  --output concat-nexus \
  --nexus
```

## 18. 版本状态

当前 workspace 版本为 `0.1.0`，Rust edition 为 2021。建议在正式数据分析中记录：

```bash
phyluce --version
phyluce config inspect
```

并保存每一步命令的 `--log-path` 日志，以便复现。

