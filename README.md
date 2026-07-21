# phyluce Rust CLI

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

- `probe easy-stampy`：用 [probebwa](https://github.com/GUIBA-EX/probebwa)
  （stampy 算法的独立 Rust 复刻）替代教程里手动调用的 `stampy.py`，一条命令
  跑完建索引、建哈希表、比对三步；已有索引时自动跳过重建，`--bam` 直接产出
  BAM，不用再手动 `samtools view`。
- `merge-multiple-gzip-files --trimmed` 和 `rename-tree-leaves --reroot`：
  原版这两个选项存在但功能缺失，这个版本补齐了。

**性能上的差异**（不改变输出，只是更快）：

移植过程中做了一批实测过的性能优化——`fast_extract` 手写扫描把默认 `--regex`
下的 contig/probe 匹配提速约 2.7x；`concatenate` 等命令的 taxon 匹配从
O(n²) 线性扫描改成哈希查找；三处逐行 SQLite INSERT 补上显式事务（单条最多
~700x）；并发调度从手撸线程池换成 `rayon`；`[profile.release]` 开启 LTO、
`phyluce-cli` 换用 `mimalloc`。这些改动都不影响任何命令的输出内容，纯粹是
内部实现，细节和实测数字见提交历史与代码注释。

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

本仓库（Rust 移植版本）的著作权归 GUIBA-EX 所有；原版 Python `phyluce` 的著作权
归 Brant C. Faircloth 所有——具体条款见 [LICENSE](LICENSE)（BSD-3-Clause，
保留原作者版权声明）。

如在论文中使用本软件，请引用本仓库；对应的 bioRxiv 预印本正在准备中，发布后
会更新此处的引用信息。机器可读的引用元数据见 [CITATION.cff](CITATION.cff)。
