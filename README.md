# phyluce Rust CLI

[![CI](https://github.com/GUIBA-EX/phyluce_rust/actions/workflows/ci.yml/badge.svg)](https://github.com/GUIBA-EX/phyluce_rust/actions/workflows/ci.yml)

[phyluce](https://github.com/faircloth-lab/phyluce)（UCE 系统发育基因组学工具包）
的 Rust 移植版本：同一套命令、同一套旧脚本名，单一静态二进制，不需要 Python/conda
环境。

## 这是什么

原版 phyluce 是 74 个独立 Python 脚本，各自依赖 Biopython/dendropy 等一整套
Python 环境。本仓库把这 74 个脚本的功能全部移植成一个 Rust 二进制
`phyluce`，按 `<domain> <command>` 分组调用：

```bash
phyluce align convert-degen-bases --alignments in --output out
```

旧脚本名依然可用——通过旧名字的 symlink 或复制后的可执行文件调用会自动映射到
新命令，例如 `phyluce_align_convert_degen_bases` 等价于
`phyluce align convert-degen-bases`。完整映射表见
[rust-command-compatibility.md](rust-command-compatibility.md)。

外部程序（MAFFT、LASTZ、SPAdes 等）仍按需调用，路径通过 `phyluce.conf` 配置；
原始 reads 的接头/质量清理应在进入本 CLI 前完成，这一点跟原版一致。

## 跟原版的区别

**行为上不一样的地方**（会影响输出或使用方式）：

- `match-contigs-to-barcodes` 不做 BOLD 数据库网络查询；改成本地 LASTZ slicing，
  传 `--no-bold` 使用。
- bootstrap replicate 用纯文本格式，不是原版的 Python `pickle`；同一个流程里
  两个实现的中间文件不能混用。
- 少数原版脚本在现代环境下本来就跑不起来（比如 Python 2 遗留代码、被移除的
  Biopython API）。这些命令按"脚本本来想干什么"实现，而不是复现它们的运行时
  报错。
- 涉及随机/并列选择的地方（tie-breaking、抽样）改成确定性规则或显式种子，
  避免原版依赖不可控随机状态导致结果不可复现。
- 部分历史遗留的 alignment 输出格式、genetree 树文件格式尚未移植，对应选项会
  明确报错，不会静默改变行为。

**原版没有、这个版本新增的命令**：

- `probe easy-stampy`：用 [probebwa](https://github.com/GUIBA-EX/probebwa)替代教程里手动调用的 `stampy.py`，一条命令
  跑完建索引、建哈希表、比对三步；已有索引时自动跳过重建，`--bam` 直接产出
  BAM，不用再手动 `samtools view`。
- `merge-multiple-gzip-files --trimmed` 和 `rename-tree-leaves --reroot`：
  原版这两个选项存在但功能缺失，这个版本补齐了。

性能优化（更快，但不改变任何命令的输出）单开一节，见下方
[性能优化](#性能优化)。

## 目录结构

- `crates/phyluce-cli`：`phyluce` 可执行文件和命令行入口。
- `crates/phyluce-align`：比对文件解析、写出、修剪、拼接和位点统计。
- `crates/phyluce-assembly`：assembly 与 match-count 相关共享逻辑。
- `crates/phyluce-io`：FASTA/FASTQ、LASTZ、2bit，以及 SQL 辅助函数。
- `crates/phyluce-config`：`phyluce.conf` 查找与 `$CONDA`/`$WORKFLOWS` 展开。
- `crates/phyluce-external`：外部命令执行封装。
- `crates/phyluce-genetrees`：gene tree 相关辅助逻辑。
- `tests/compat`：兼容性和回归测试脚本。

文档入口见 [docs/index.md](docs/index.md)，完整中文手册见
[docs/user-manual-zh.md](docs/user-manual-zh.md)（构建配置、端到端 UCE 流程、
各命令域示例、故障排查）；已知差异的完整列表见
[rust-command-compatibility.md](rust-command-compatibility.md)。

## 性能优化

原则是"先跑 benchmark，能证明有收益再改"——下表只列有实测数字支撑的优化，
且全部不改变任何命令的输出内容，纯粹是内部实现；细节见提交历史和对应代码里
的 benchmark 注释。

| 优化点 | 方法 | 实测提升 |
| --- | --- | --- |
| contig/probe 名称匹配（默认 `--regex`） | 手写扫描替代通用正则引擎（`fast_extract`），配差分模糊测试防止跟正则语义跑偏 | ~2.7x |
| `concatenate` 等命令的 taxon 匹配 | O(n²) 线性扫描 → 哈希查找 | 消除随 taxon 数二次方增长（`concatenate` 等 3 处） |
| 逐行 SQLite INSERT（3 处命令） | 补上显式事务，避免每行触发一次 autocommit | 最多 ~700x |
| 并发任务调度 | `rayon` | ~1.8x–4.4x |
| FASTA 解析 | 消除序列行的重复扫描（原来每行最多扫 3 遍） | ~1.3x |
| FASTQ 长度/计数 | 逐行 `BufRead::lines()`（每行分配 `String` + UTF-8 校验）→ 字节级 `read_until` 读入复用缓冲区，只在真正需要时才转 UTF-8 | 长度提取 ~1.2x，纯计数 ~3.3x |
| 2bit 解码 | 逐 base 除法+移位+查表 → 逐字节查 256 项表（一次查表拿到 4 个 base，`extend_from_slice` 写入） | ~2.8x–3x（~1.2 → ~3.3-3.6 Gbases/sec） |
| 编译配置 | `[profile.release]` 开启 LTO + `codegen-units=1` | - |
| 内存分配器 | `phyluce-cli` 换用 `mimalloc` | - |

**验证过但没有采用的优化**（同样先测再决定，负结果也如实记录）：

- `ahash` 替代标准库 SipHash：已实测，收益边际，仅在已经因为改用哈希表而
  受益的地方顺手用上，没有为此单独改动。
- SIMD 解码 2bit：最初只测了朴素的逐 base 循环（~1.2 Gbases/sec），据此认为
  "已经够快+`std::simd` 是 nightly-only"就没再深入——这个结论下得太早，见上表，
  一个不需要 SIMD/unsafe/nightly 的查表方案就有 ~3x 提升，反而佐证了"先把
  简单方案的账算清楚，再谈需不需要 SIMD"。查表方案封顶后（~3.5 Gbases/sec，
  解码一整条人类染色体不到一秒）已经远超本项目实际需求，SIMD 在这基础上还能
  再快多少没有验证，暂不追加。
- `mmap` 读取 2bit 文件：实测评估后判断收益不足以抵消复杂度。
- 引入 [rust-bio](https://github.com/rust-bio/rust-bio)（crates.io 上的
  `bio`）替代自己的 FASTA/FASTQ reader：单独测过，它确实更快（同一份数据
  上 FASTA 快 ~2x、FASTQ 快 ~1.8x），但会带来 83 个传递依赖（`ndarray`、
  `nalgebra`、`statrs`、`petgraph` 等一整套统计/图/线性代数库），而
  `phyluce-io` 是全项目每个命令都依赖的基础 crate；真实 UCE 流程里外部
  工具（LASTZ/MAFFT/SPAdes 等）的子进程耗时才是大头，不是本身的文件解析。
  没有加这个依赖，改成照着差距手动优化了上面那两行。

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

兼容性测试（只用仓库内 fixture，可在独立克隆中运行）：

```bash
python3 tests/compat/run_fixtures.py target/debug/phyluce
```

要跑包含原版 Python 和外部工具的完整对照，设置
`PHYLUCE_PYTHON_REPO=/path/to/phyluce` 后运行
`python3 tests/compat/run_all.py target/debug/phyluce`。

## 版权与引用

本仓库（Rust 移植版本）的著作权归 GUIBA-EX 所有，采用 **GPLv3
（GPL-3.0-or-later）** 许可；原版 Python `phyluce` 的著作权归 Brant C.
Faircloth 所有（BSD-3-Clause）——原始版权声明按其许可证要求保留在
[LICENSE](LICENSE) 文件中。

如在论文中使用本软件，请引用本仓库；对应的 bioRxiv 预印本正在准备中，发布后
会更新此处的引用信息。机器可读的引用元数据见 [CITATION.cff](CITATION.cff)。
