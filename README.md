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

## 快速开始

```bash
cargo build -p phyluce-cli --release
target/release/phyluce --version
target/release/phyluce --help
target/release/phyluce config inspect
```

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
