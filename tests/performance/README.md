# 比对与拼接性能测试

`benchmark_alignment.py` 对比优化前后 release 二进制，使用真实 MAFFT，检查全部输出
文件的 SHA-256 和标准输出一致。每个场景运行三次，交替版本顺序，不清理系统缓存。
矩阵输入为 1000 loci × 200 taxa × 500 bp、80% occupancy；最终矩阵 100 MB。
MAFFT 使用 4 个 worker：均匀输入为 20 loci × 8 taxa × 400 bp；不均匀输入将最后
4 个 locus 增为 24 taxa × 1800 bp。序列由固定随机种子生成，包含约 6% 随机替换。

2026-10-19，Linux x86_64，MAFFT 7.505，中位数结果：

| 场景 | 优化前 | 优化后默认 | 耗时变化 |
|---|---:|---:|---:|
| NEXUS 拼接 | 190.6 ms | 95.2 ms | −50.1% |
| PHYLIP 拼接 | 177.9 ms | 89.0 ms | −50.0% |
| MAFFT 均匀、不修剪 | 347.9 ms | 340.6 ms | −2.1% |
| MAFFT 均匀、修剪 | 351.3 ms | 352.3 ms | +0.3% |
| MAFFT 不均匀、不修剪 | 586.3 ms | 525.4 ms | −10.4% |
| MAFFT 不均匀、修剪 | 580.5 ms | 525.3 ms | −9.5% |

拼接峰值 RSS 从约 24.7 MiB 增至 212.7 MiB，以内存换取速度。候选版本强制磁盘
路径的 NEXUS/PHYLIP 耗时为 197.4/190.7 ms，比基线慢约 3.6%/7.2%，RSS 约
18.5 MiB。100 MiB 预算的矩阵驻内存、输入重新读取路径耗时为 110.2/106.6 ms，
RSS 约 162.7 MiB；释放缓存后分配器可能保留内存，预算不等于 RSS 上限。

这些是小规模合成数据的阶段测试，均匀任务的差异很小，不代表可靠加速；没有声称
完整工作流或所有真实数据均可获得相同收益。GNU time 的 RSS 不代表多个 MAFFT
子进程的同时驻留内存之和。原始各次耗时、CPU 时间、RSS 和二进制指纹保存在
[结果 JSON](results/2026-10-19-alignment.json)。

复现（需要 GNU time 和 MAFFT）：

```bash
python3 tests/performance/benchmark_alignment.py \
  --baseline /path/to/before/phyluce \
  --candidate target/release/phyluce \
  --mafft /path/to/mafft \
  --output /tmp/alignment-benchmark.json
```

脚本会生成临时配置选择 MAFFT。若使用解包但未安装的 MAFFT，需额外设置
`MAFFT_BINARIES` 指向其辅助程序目录。所有输入与输出置于自动清理的临时目录。

此前的提取/匹配测试由 `benchmark_pipeline.py` 复现，原始结果位于
`results/2026-10-19.json`。

## 后续内存优化：直接输出缓存片段

`results/2026-10-19-concat-memory.json` 的基线为上文已经优化过的版本，候选版本
不再同时分配缓存和完整矩阵。相同的 100 MB 矩阵输入、五次交替运行中位数：

| 输出 | 上一版耗时 | 本版耗时 | 上一版 RSS | 本版 RSS |
|---|---:|---:|---:|---:|
| NEXUS | 89.1 ms | 66.9 ms | 212.8 MiB | 120.5 MiB |
| PHYLIP | 87.8 ms | 66.1 ms | 212.7 MiB | 120.5 MiB |

默认路径峰值 RSS 减少约 43%，耗时减少约 25%。100 MiB 预算现在也能保留输入
缓存，耗时约 66–68 ms、RSS 约 121 MiB（预算不是 RSS 硬上限）。强制磁盘路径
耗时约 184–187 ms、RSS 约 18.6 MiB。所有版本、预算及格式的输出逐字节一致。
上文历史结果中的 100 MiB 路径描述仅适用于上一版。

复现使用同一脚本，将 `--baseline` 指向上一版二进制，并增加
`--concat-only --repeats 5`。合成数据和热缓存的限制仍适用。

## 解析后立即打包缓存、合并序列分配

基线为上一节的直接输出缓存版本。每个 locus 解析后立即用整数 taxon 编号替代
重复名称，并将序列合并到一个连续缓冲区；输出片段索引按实际条数预分配。
预算估算包括序列、taxon 编号及后续片段索引，仍不是 RSS 硬上限。

测试使用长 taxon 名称、不同 taxon/locus 的序列内容及 N、缺口、缺失字符。
常规数据为 1000 loci × 200 taxa × 500 bp、80% occupancy；短 locus 数据为
5000 × 200 × 40 bp、80% occupancy；高缺失率数据为 1000 × 200 × 500 bp、
10% occupancy。热缓存、五次交替运行中位数如下（默认 256 MiB 预算）：

| 场景 | 格式 | 基线耗时 | 本版耗时 | 基线 RSS | 本版 RSS |
|---|---|---:|---:|---:|---:|
| 常规 | NEXUS | 74.3 ms | 71.2 ms | 120.9 MiB | 108.6 MiB |
| 常规 | PHYLIP | 73.9 ms | 70.8 ms | 120.9 MiB | 108.6 MiB |
| 大量短 locus | NEXUS | 163.3 ms | 151.9 ms | 216.8 MiB | 104.6 MiB |
| 大量短 locus | PHYLIP | 162.9 ms | 153.7 ms | 216.7 MiB | 104.5 MiB |
| 高缺失率 | NEXUS | 30.8 ms | 32.1 ms | 26.9 MiB | 26.5 MiB |
| 高缺失率 | PHYLIP | 30.8 ms | 31.6 ms | 26.8 MiB | 26.6 MiB |

短 locus 的峰值内存约减半，常规数据降低约 10%。高缺失率场景收益很小，
默认路径耗时增加约 0.9–1.4 ms，不能声称所有场景都加速。32 MiB 预算下，常规
和短 locus 场景仍会回退，耗时没有稳定改善；本次未改变超预算后重新解析的策略。

共 150 次运行覆盖 NEXUS/PHYLIP、基线及本版默认路径、双方 32 MiB 预算和本版
强制磁盘路径；输出文件及标准输出全部一致。已有跨片段、缺失 taxon、磁盘回退
测试以及上游 NEXUS/PHYLIP/FASTA 输入兼容性测试通过。原始结果见
[packed-cache JSON](results/2026-10-19-packed-cache.json)。

```bash
python3 tests/performance/benchmark_packed_cache.py \
  --baseline /path/to/previous/phyluce \
  --candidate target/release/phyluce \
  --repeats 5 \
  --output /tmp/packed-cache.json
```

此测试只依赖 GNU time，不需要 MAFFT；合成数据的阶段测试不代表完整工作流性能。

## 超预算回退：单次解析与二进制快照

基线为上一节的紧凑缓存版本。本版在超预算时转存已有缓存，并将后续解析结果
追加为二进制记录（taxon 编号、矩阵起点、长度、序列）；构建矩阵时顺序读取快照，
不再读取或解析原始输入。原先用于两遍输入一致性检查的逐 locus taxon 集合也已
移除。原始输入仍应在读取期间保持稳定；输出来自各文件首次成功解析的内容。

快照和矩阵均使用自动删除的临时文件。快照截断会报错，成功及失败路径都会清理；
单元测试还验证了原始文件删除后仍能使用快照完成拼接，以及缓存中途超预算。
快照占用实际序列字节数加每条序列 24 字节，构建磁盘矩阵时两者同时占用磁盘。
预算仍不包含解析器、元数据和内存分配器保留空间，不是进程 RSS 硬上限。

FASTA 输入沿用前一节三类数据，五次交替运行中位数（下表为 NEXUS 输出）：

| 场景 | 预算 | 基线耗时 | 本版耗时 | 基线 RSS | 本版 RSS |
|---|---:|---:|---:|---:|---:|
| 常规 | 0 MiB | 195.9 ms | 188.5 ms | 16.8 MiB | 14.7 MiB |
| 常规 | 32 MiB | 200.7 ms | 197.3 ms | 48.9 MiB | 46.6 MiB |
| 常规 | 48 MiB | 202.6 ms | 195.2 ms | 65.1 MiB | 62.7 MiB |
| 大量短 locus | 0 MiB | 601.3 ms | 474.1 ms | 38.9 MiB | 22.8 MiB |
| 大量短 locus | 32 MiB | 611.8 ms | 478.8 ms | 54.9 MiB | 44.7 MiB |
| 大量短 locus | 48 MiB | 279.7 ms | 162.6 ms | 110.8 MiB | 90.7 MiB |
| 高缺失率 | 0 MiB | 63.7 ms | 62.5 ms | 14.8 MiB | 14.7 MiB |
| 高缺失率 | 32 MiB | 31.1 ms | 30.0 ms | 26.9 MiB | 24.7 MiB |
| 高缺失率 | 48 MiB | 30.9 ms | 30.1 ms | 26.8 MiB | 24.5 MiB |

常规和短 locus 的 32 MiB 路径使用磁盘矩阵；短 locus 的 48 MiB 路径使用
内存矩阵。高缺失率的 32/48 MiB 路径仍能容纳缓存，不触发回退。短 locus 回退
收益明显，常规和高缺失率数据收益较小。热缓存结果不代表慢磁盘或完整工作流。

FASTA 共 240 次运行，覆盖 NEXUS/PHYLIP 输出、双方默认/0/32/48 MiB 预算；
全部输出文件及标准输出一致。原始结果见
[spool JSON](results/2026-10-19-spool.json)。

```bash
python3 tests/performance/benchmark_packed_cache.py \
  --baseline /path/to/packed-cache/phyluce \
  --candidate target/release/phyluce --repeats 5 \
  --output /tmp/spool.json
```

补充 NEXUS 输入测试：5000 个短 locus，三次交替运行，双方四种预算、两种输出
共 48 次运行全部一致。32 MiB 预算、NEXUS 输出耗时从 809.5 ms
降至 572.8 ms，RSS 从 56.8 MiB 降至 44.5 MiB。
原始结果见 [NEXUS 输入 JSON](results/2026-10-19-spool-nexus.json)。
复现时在上述命令增加 `--input-format nexus --cases short-loci --repeats 3`。

## 直接从磁盘快照输出，取消完整磁盘矩阵

基线为上一节单次解析快照版本。本版磁盘路径只顺序追加实际序列，内存中保留
每条序列的起点、长度和文件偏移。输出时按 taxon 读取所需片段，缺失位置即时
补 `?`，不再初始化或随机写入完整磁盘矩阵。Unix 使用定位读取，避免每条片段
先 seek 再 read。矩阵能放进预算时仍保留内存路径，按文件偏移排序后顺序读取
快照；默认缓存足够时的路径不变。

五次交替运行中位数；下表为 FASTA 输入、NEXUS 输出、强制磁盘模式：

| 场景 | 基线耗时 | 本版耗时 | 临时磁盘峰值（基线→本版） | RSS（基线→本版） |
|---|---:|---:|---:|---:|
| 常规 | 165.4 ms | 117.7 ms | 183.84 → 80.00 MB | 15.0 → 18.8 MiB |
| 大量短 locus | 345.4 ms | 246.9 ms | 91.20 → 32.00 MB | 23.0 → 53.1 MiB |
| 高缺失率 | 60.0 ms | 36.6 ms | 110.48 → 10.00 MB | 14.9 → 14.8 MiB |

代价是索引驻留内存。短 locus 的 48 MiB 内存矩阵路径耗时从 159.6 ms
增至 169.5 ms，RSS 从 90.9 MiB 增至 131.0 MiB；其临时
占用仍从 51.2 MB 降至 32 MB。本次侧重磁盘写入和占用，并非所有路径都加速。
索引及排序引用不包含在矩阵预算中，预算不是 RSS 硬上限。

测试增加隔离的 TMPDIR，以约 1 ms 间隔采样临时文件逻辑长度和实际分配空间，
并断言每次运行后临时文件全部清理。表中为采样到的逻辑长度峰值，不包括最终
输出；采样可能漏掉极短瞬时峰值。观测值与本组数据的预期快照长度一致。实际
分配空间原始值见 JSON。监测开销同时作用于两版。仍是热缓存阶段测试，机械
硬盘或冷缓存下的随机读取性能尚未验证。

FASTA 共 240 次运行，覆盖三类输入、两种输出及双方四种预算，所有输出文件和
标准输出完全一致。跨 taxon/locus 的部分读取、负偏移、快照截断和清理测试均
通过。结果见 [磁盘输出 JSON](results/2026-10-19-disk-output.json)。

```bash
python3 tests/performance/benchmark_packed_cache.py \
  --baseline /path/to/spool-version/phyluce \
  --candidate target/release/phyluce --repeats 5 --track-temp-disk \
  --output /tmp/disk-output.json
```

补充 NEXUS 输入：5000 个短 locus、三次交替运行，共 48 次输出对照全部一致。
强制磁盘、NEXUS 输出耗时从 436.1 ms 降至 338.9 ms，
临时占用同样从 91.2 MB 降至 32 MB。完整预算/格式/RSS 数据见
[NEXUS 输入结果](results/2026-10-19-disk-output-nexus.json)。
复现时增加 `--input-format nexus --cases short-loci --repeats 3`。
