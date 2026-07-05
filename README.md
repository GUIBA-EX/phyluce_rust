# phyluce Rust CLI

本目录包含 `phyluce` Rust CLI 的 workspace。代码按功能拆分为多个 crate，
由 `phyluce-cli` 统一暴露命令入口。

## 目录结构

- `crates/phyluce-cli`：`phyluce` 可执行文件和命令行入口。
- `crates/phyluce-align`：比对文件解析、写出、修剪、拼接和位点统计。
- `crates/phyluce-assembly`：assembly 与 match-count 相关共享逻辑。
- `crates/phyluce-io`：FASTA/FASTQ、LASTZ、2bit，以及 SQL 辅助函数。
- `crates/phyluce-config`：`phyluce.conf` 查找与 `$CONDA`/`$WORKFLOWS` 展开。
- `crates/phyluce-external`：外部命令执行封装。
- `crates/phyluce-genetrees`：gene tree 相关辅助逻辑。
- `tests/compat`：兼容性和回归测试脚本。

## 用户手册

中文用户手册见 [docs/user-manual-zh.md](docs/user-manual-zh.md)。手册按原版
`phyluce_*` 脚本的使用方式组织，包含安装配置、UCE 主流程、各命令域示例、
旧脚本名兼容、已知差异和故障排查。

## CLI 形式

Rust 版本使用一个分组式 CLI：

```text
phyluce <domain> <command> [options]
```

示例：

```text
phyluce align convert-degen-bases --alignments in --output out
```

当前这批 align/assembly 命令也支持旧脚本名兼容：如果通过旧脚本名的
symlink 或复制后的可执行文件调用，会自动映射到新的分组命令。例如
`phyluce_align_convert_degen_bases` 会映射为：

```text
phyluce align convert-degen-bases
```

完整旧命令名映射和日志行为见
[rust-command-compatibility.md](rust-command-compatibility.md)。

## 主要改动

- 新增 16 个 align/assembly 命令移植，包括 degenerate-base 转换、
  alignment explode、taxon FASTA 提取、PAML 分区格式化、matrix estimates、
  min-taxa 过滤、taxon-locus 计数、按配置移动 alignment、随机抽样拼接、
  RAxML reduction、移除 locus name、problem screening、smilogram、probe
  duplicate screening、barcode extraction 和 barcode matching。
- 修复格式识别，避免将 `phylip`、`phylip-relaxed`、`clustal`、`emboss`、
  `stockholm` 静默当作 FASTA 处理。
- 新增 `phyluce-io::sql`，集中处理 SQL 标识符转义，并替换相关动态表名、
  列名拼接。
- 新增旧命令名兼容层，方便现有命令调用方式逐步迁移。
- 新增基于 `tracing` 的可选文件日志：全局参数 `--log-path` 和
  `--verbosity`。默认不写日志，不改变 stdout/stderr。
- 扩展兼容性测试，优先使用已有 fixture，随机、外部工具或历史兼容问题路径
  保留合成 smoke test。

## 已知差异

- `match-contigs-to-barcodes` 不执行 BOLD 网络查询；需要本地 LASTZ slicing
  时请传入 `--no-bold`。
- 依赖外部工具的命令仍需要正确配置 MAFFT、LASTZ、RAxML 等路径。
- 少数历史脚本存在运行时兼容性问题；Rust 版本按预期行为实现，而不是复现
  运行时失败。

## 开发检查

在 `rust/` 目录下运行：

```text
cargo check -p phyluce-cli
cargo test -p phyluce-io -p phyluce-assembly -p phyluce-cli
cargo clippy -p phyluce-cli --all-targets -- -D warnings
cargo build -p phyluce-cli
```

兼容性测试：

```text
python3 tests/compat/compare_align_assembly_batch2.py target/debug/phyluce
python3 tests/compat/run_all.py target/debug/phyluce
```

`run_all.py` 会调用 MAFFT 等外部工具；在受限 sandbox 中可能需要在非 sandbox
环境运行。
