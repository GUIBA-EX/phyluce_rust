[![CI](https://github.com/GUIBA-EX/phyluce_rust/actions/workflows/ci.yml/badge.svg)](https://github.com/GUIBA-EX/phyluce_rust/actions/workflows/ci.yml)

# phyluce_rust

<p align="center">
  <img src="docs/assets/logo.png" alt="phyluce_rust logo" width="240">
</p>

[phyluce](https://github.com/faircloth-lab/phyluce)（UCE 系统发育基因组学工具包）的
Rust 移植版本：命令集与旧脚本名称保持一致，编译为单一静态二进制文件，无需
Python 或 conda 环境。

## 简介

原版 phyluce 由 74 个独立的 Python 脚本组成，各自依赖 Biopython、DendroPy 等
完整 Python 环境。本仓库将其全部功能移植为单一 Rust 二进制 `phyluce`，按
`<domain> <command>` 分组调用：

```bash
phyluce align convert-degen-bases --alignments in --output out
```

旧脚本名称仍可使用：以旧名称建立符号链接或复制可执行文件，调用时会自动映射至
新命令，例如 `phyluce_align_convert_degen_bases` 等价于
`phyluce align convert-degen-bases`。完整映射表见
[rust-command-compatibility.md](rust-command-compatibility.md)。

外部程序（MAFFT、LASTZ、SPAdes 等）仍按需调用，路径由 `phyluce.conf` 配置；
原始 reads 的接头与质量清理须在进入本 CLI 前完成，此点与原版一致。

## 与原版的差异

**行为差异**（影响输出或使用方式）：

- `match-contigs-to-barcodes` 不执行 BOLD 数据库网络查询，改为本地 LASTZ
  切片，须传入 `--no-bold`。
- bootstrap replicate 采用纯文本格式，而非原版的 Python `pickle`；两种实现
  产生的中间文件不可混用。
- `prep-uce-align-files-for-ncbi`（对应原版
  `phyluce_ncbi_prep_uce_align_files_for_ncbi`）在现代 Biopython 环境下因
  `Bio.Alphabet` 模块已被移除而无法导入运行；本版本按其预期行为实现，而非
  复现该运行时错误。
- 涉及随机或并列选择的逻辑（如 tie-breaking、抽样）改为确定性规则或显式随机
  种子，避免原版因不可控随机状态导致结果不可复现。
- 部分历史遗留的 alignment 输出格式与 genetree 树文件格式尚未移植，相关选项
  会明确报错，不会静默改变行为。

**新增命令**（原版没有）：

- [`probe get-full-probes`](#生成-geneminer2-参考序列)：从原版 probe FASTA
  生成 GeneMiner2 逐-locus 参考目录。
- `probe easy-stampy`：以 [probebwa](https://github.com/GUIBA-EX/probebwa)
  替代教程中手动调用的 `stampy.py`，一条命令完成建索引、建哈希表、比对三
  步；索引已存在时自动跳过重建，`--bam` 直接输出 BAM，无需再手动调用
  `samtools view`。
- `merge-multiple-gzip-files --trimmed` 与 `rename-tree-leaves --reroot`：
  原版存在对应选项但功能缺失，本版本已补全。

性能优化详见下文[性能优化](#性能优化)一节。

## 目录结构

- `crates/phyluce-cli`：`phyluce` 可执行文件与命令行入口。
- `crates/phyluce-align`：比对文件的解析、写出、修剪、拼接与位点统计。
- `crates/phyluce-assembly`：assembly 与 match-count 相关共享逻辑。
- `crates/phyluce-io`：FASTA/FASTQ、LASTZ、2bit 及 SQL 辅助函数。
- `crates/phyluce-config`：`phyluce.conf` 查找与 `$CONDA`/`$WORKFLOWS` 展开。
- `crates/phyluce-external`：外部命令执行封装。
- `crates/phyluce-genetrees`：gene tree 相关辅助逻辑。
- `tests/compat`：兼容性与回归测试脚本。

文档入口见 [docs/index.md](docs/index.md)，完整中文手册见
[docs/user-manual-zh.md](docs/user-manual-zh.md)（构建配置、端到端 UCE
流程、各命令域示例、故障排查）；已知差异的完整列表见
[rust-command-compatibility.md](rust-command-compatibility.md)。

## 性能优化

指导原则：仅在有实测数据支持时引入优化，且优化不改变任何命令的输出内容。
下表所列优化均有具体数字支撑，方法细节见对应代码中的 benchmark 注释与提交
历史。

| 优化点 | 方法 | 实测提升 |
| --- | --- | --- |
| contig/probe 名称匹配（默认 `--regex`） | 手写扫描替代通用正则引擎（`fast_extract`） | ~2.7x |
| `concatenate` 等命令的 taxon 匹配 | O(n²) 线性扫描改为哈希查找 | 消除二次方增长 |
| 逐行 SQLite INSERT（3 处命令） | 补充显式事务 | 最多 ~700x |
| 并发任务调度 | 改用 `rayon` | ~1.8x–4.4x |
| FASTA/FASTQ 解析 | 消除逐行重复扫描；FASTQ 改用字节级 `read_until` | FASTA ~1.3x，FASTQ 长度提取 ~1.2x、计数 ~3.3x |
| 2bit 解码 | 逐 base 查表改为逐字节查 256 项表 | ~3x（约 3.5 Gbases/秒） |
| 编译配置与内存分配器 | LTO + `codegen-units=1`，改用 `mimalloc` | - |

### 流程级优化与基准

`assembly match-contigs-to-probes` 新增 `--cores N`（默认 1）进行有界样本
并行，新增 `--resume` 校验并复用已完成的 LASTZ 结果。复用以输入、可执行
文件、参数及输出内容指纹为依据；每次重新生成 SQLite 和报告。
`get-fastas-from-match-counts` 和 FASTA 长度统计改为流式读取；提取输出使用
缓冲写入。匹配数据库更新在每个样本事务内复用预编译 SQL 语句。

2026-10-19 本地 release 基准（5 次运行中位数）：

| 阶段 | 优化前 | 优化后 |
| --- | --- | --- |
| 提取峰值 RSS（250,000 条 × 1,000 bp，约 258 MB FASTA） | 300.7 MiB | 22.5 MiB |
| 同一提取任务耗时 | 84.6 ms | 66.1 ms |
| 匹配结果处理（64 样本 × 5,000 loci，单 worker） | 2.22 s | 1.46 s |

这是合成数据、主要命中缓存的阶段测量，没有执行 LASTZ/MAFFT/组装器；结果
处理的四 worker 测量为 1.51 s，未优于单 worker，默认并行度因此保持 1。
输出 FASTA、SQLite 内容和 CSV 已核对一致。基准脚本见
[tests/performance/benchmark_pipeline.py](tests/performance/benchmark_pipeline.py)，
原始测量、二进制校验和及环境信息见
[2026-10-19.json](tests/performance/results/2026-10-19.json)。

## 快速开始

```bash
cargo build -p phyluce-cli --release
target/release/phyluce --version
target/release/phyluce --help
target/release/phyluce config inspect
```

## 生成 GeneMiner2 参考序列

```bash
phyluce probe get-full-probes \
  --input probes.fasta \
  --output ref_uce
```

输入须为原版 probe FASTA，header 中包含 local/global 坐标，序列仅含 A/C/G/T。
命令按坐标拼接重叠 probe，每组保留最长有效连续序列，并在输出目录中为每个
locus 写一个 FASTA 文件。输出目录须不存在或为空。不使用 MUSCLE 或 MAFFT。

## 比对与拼接性能

MAFFT 多进程任务优先调度估算工作量较大的 locus，保留原有比对参数和输出顺序。
矩阵拼接从缓存的序列片段直接输出，避免缓存与完整矩阵同时驻留；默认提供
256 MiB 的矩阵或输入缓存预算；`--memory-limit-mb 0` 可强制
使用低内存磁盘路径。超出缓存预算后转存二进制临时快照，避免重新解析原始输入；
磁盘路径直接读取快照输出，不再构建完整磁盘矩阵。临时文件只保存实际序列，
片段索引驻留内存；预算不是进程总内存上限。

实测结果及复现方式见 [比对与拼接性能测试](tests/performance/README.md)。

## 开发检查

```bash
cargo check -p phyluce-cli
cargo test -p phyluce-io -p phyluce-assembly -p phyluce-cli
cargo clippy -p phyluce-cli --all-targets -- -D warnings
cargo build -p phyluce-cli
```

兼容性测试（仅使用仓库内置 fixture，可在独立克隆中运行）：

```bash
python3 tests/compat/run_fixtures.py target/debug/phyluce
```

如需包含原版 Python 及外部工具的完整对照测试，设置
`PHYLUCE_PYTHON_REPO=/path/to/phyluce` 后运行
`python3 tests/compat/run_all.py target/debug/phyluce`。

## 版权与引用

本仓库（Rust 移植版本）的著作权归 GUIBA-EX 所有，采用 **GPLv3
（GPL-3.0-or-later）** 许可；原版 Python `phyluce` 的著作权归 Brant C.
Faircloth 所有（BSD-3-Clause）——原始版权声明按其许可证要求保留于
[LICENSE](LICENSE) 文件中。

若在论文中使用本软件，请引用本仓库；对应的 bioRxiv 预印本正在准备中，发布
后将更新此处的引用信息。机器可读的引用元数据见 [CITATION.cff](CITATION.cff)。
