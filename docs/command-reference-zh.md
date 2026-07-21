# 命令速查表

按命令域列出全部命令：新形式（`phyluce <domain> <command>`）、对应的旧脚本名，
以及一句话说明。说明文字参照
[原版 PHYLUCE 的 List of Programs](https://phyluce.readthedocs.io/en/latest/daily-use/list-of-programs.html)
整理，按 Rust 版当前行为校对；某条命令的完整参数请用
`phyluce <domain> <command> --help` 查看，用法示例见
[用户手册](user-manual-zh.md)对应章节，已知差异见
[rust-command-compatibility.md](../rust-command-compatibility.md)。

## align（27 个）

| 新命令 | 旧脚本名 | 说明 |
| --- | --- | --- |
| `align add-missing-data-designators` | `phyluce_align_add_missing_data_designators` | 给缺少某些 taxa 的 alignment 补上缺失标记；用 phyluce 的拼接工具时通常不需要单独跑这个。 |
| `align concatenate-alignments` | `phyluce_align_concatenate_alignments` | 把一批 alignment 拼接成一个 NEXUS 或 PHYLIP 文件，附带 charset 信息。 |
| `align convert-degen-bases` | `phyluce_align_convert_degen_bases` | 把 alignment 里的 IUPAC 简并碱基码统一转成 `N`。 |
| `align convert-one-align-to-another` | `phyluce_align_convert_one_align_to_another` | 在 FASTA/NEXUS/PHYLIP/PHYLIP-relaxed/CLUSTAL/EMBOSS/Stockholm 之间自由转换格式。 |
| `align explode-alignments` | `phyluce_align_explode_alignments` | 把一个目录的 alignment "炸开"成按 taxon 或按 locus 拆分的独立文件。 |
| `align extract-taxa-from-alignments` | `phyluce_align_extract_taxa_from_alignments` | 按指定 taxa 列表保留或剔除，生成新的 alignment 目录。 |
| `align extract-taxon-fasta-from-alignments` | `phyluce_align_extract_taxon_fasta_from_alignments` | 从一批 alignment 里抽出某个 taxon 的数据，输出为 FASTA。 |
| `align filter-alignments` | `phyluce_align_filter_alignments` | 按 taxa 数量或长度筛选 alignment，把不合格的排除到新目录之外。 |
| `align format-concatenated-phylip-for-paml` | `phyluce_align_format_concatenated_phylip_for_paml` | 把拼接好的 PHYLIP alignment 转成 PAML 需要的内部格式。 |
| `align get-align-summary-data` | `phyluce_align_get_align_summary_data` | 快速汇总一批 alignment 的统计信息。 |
| `align get-gblocks-trimmed-alignments-from-untrimmed` | `phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed` | 用 Gblocks 修剪 alignment 边缘，输出到新目录。 |
| `align get-incomplete-matrix-estimates` | `phyluce_align_get_incomplete_matrix_estimates` | 估算不同"允许缺失"阈值下矩阵里各能保留多少 taxa。 |
| `align get-informative-sites` | `phyluce_align_get_informative_sites` | 统计每个 alignment 的信息位点数，可用于按位点数进一步过滤。 |
| `align get-only-loci-with-min-taxa` | `phyluce_align_get_only_loci_with_min_taxa` | 按最小 taxa 数过滤 alignment，把满足条件的输出到新目录。 |
| `align get-ry-recoded-alignments` | `phyluce_align_get_ry_recoded_alignments` | 把 alignment 重编码成 R/Y 或 0/1 两态数据。 |
| `align get-smilogram-from-alignments` | `phyluce_align_get_smilogram_from_alignments` | 输出从 alignment 中心到边缘的位点数分布（CSV），可用来画 "smilogram" 图观察 UCE 变异模式。 |
| `align get-taxon-locus-counts-in-alignments` | `phyluce_align_get_taxon_locus_counts_in_alignments` | 统计每个 taxon 出现在多少个 alignment 里，可用于按覆盖度过滤。 |
| `align get-trimal-trimmed-alignments-from-untrimmed` | `phyluce_align_get_trimal_trimmed_alignments_from_untrimmed` | 用 trimAl 修剪 alignment，输出到新目录。 |
| `align get-trimmed-alignments-from-untrimmed` | `phyluce_align_get_trimmed_alignments_from_untrimmed` | 用 phyluce 原生的 edge-trimming 算法修剪 alignment，输出到新目录。 |
| `align move-align-by-conf-file` | `phyluce_align_move_align_by_conf_file` | 按配置文件里列出的名单，把 alignment 从一个目录复制到另一个目录，常用于筛选。 |
| `align randomly-sample-and-concatenate` | `phyluce_align_randomly_sample_and_concatenate` | 从一批 alignment 中随机抽样，把抽到的序列拼接输出。 |
| `align reduce-alignments-with-raxml` | `phyluce_align_reduce_alignments_with_raxml` | 用 RAxML 生成每个 alignment 的"精简版"（去掉缺失和低信息量的位点模式）。 |
| `align remove-empty-taxa` | `phyluce_align_remove_empty_taxa` | 去掉 alignment 里完全没有数据的 taxa。 |
| `align remove-locus-name-from-files` | `phyluce_align_remove_locus_name_from_files` | 去掉序列名里附带的 UCE locus 名称前缀/后缀，输出到新目录。 |
| `align screen-alignments-for-problems` | `phyluce_align_screen_alignments_for_problems` | 排查 alignment 里的异常情况，比如奇怪的碱基码（`X`）或连续的模糊碱基（`N`）。 |
| `align seqcap-align` | `phyluce_align_seqcap_align` | 给一个 monolithic FASTA 文件按 locus 分组做 alignment，输出到新目录（内部调用 MAFFT）。 |
| `align split-concat-nexus-to-loci` | `phyluce_align_split_concat_nexus_to_loci` | 把带 charset 信息的拼接 NEXUS 文件拆回各个 locus。 |

## assembly（13 个）

| 新命令 | 旧脚本名 | 说明 |
| --- | --- | --- |
| `assembly assemblo-abyss` | `phyluce_assembly_assemblo_abyss` | 用 ABySS 组装 fastq 数据。 |
| `assembly assemblo-spades` | `phyluce_assembly_assemblo_spades` | 用 SPAdes 组装 fastq 数据。 |
| `assembly assemblo-velvet` | `phyluce_assembly_assemblo_velvet` | 用 Velvet 组装 fastq 数据。 |
| `assembly explode-get-fastas-file` | `phyluce_assembly_explode_get_fastas_file` | 把一个 monolithic 的 UCE contig FASTA 文件按 locus 或 taxon 拆成独立文件。 |
| `assembly extract-contigs-to-barcodes` | `phyluce_assembly_extract_contigs_to_barcodes` | 把 `match-contigs-to-barcodes` 产出的日志整理成更易读的结果表格。 |
| `assembly get-bed-from-lastz` | `phyluce_assembly_get_bed_from_lastz` | 把 LASTZ 结果文件转换成 BED 格式。 |
| `assembly get-fasta-lengths` | `phyluce_assembly_get_fasta_lengths` | 汇总一个 FASTA 文件里 contig 的长度等统计信息。 |
| `assembly get-fastas-from-match-counts` | `phyluce_assembly_get_fastas_from_match_counts` | 根据 `get-match-counts` 的结果，产出一个 monolithic 的 UCE locus FASTA 文件。 |
| `assembly get-fastq-lengths` | `phyluce_assembly_get_fastq_lengths` | 汇总一批 fastq reads 的统计信息。 |
| `assembly get-match-counts` | `phyluce_assembly_get_match_counts` | 结合 `match-contigs-to-probes` 的结果和配置文件，输出各 taxon/locus 命中矩阵。 |
| `assembly match-contigs-to-barcodes` | `phyluce_assembly_match_contigs_to_barcodes` | 检查每个 taxon 的 contig 里是否含有物种条形码序列，抽取对应区域；Rust 版默认不做 BOLD 网络查询（`--no-bold`）。 |
| `assembly match-contigs-to-probes` | `phyluce_assembly_match_contigs_to_probes` | 在一批组装好的 contig 里搜索跟 UCE bait/probe 匹配的部分。 |
| `assembly screen-probes-for-dupes` | `phyluce_assembly_screen_probes_for_dupes` | 检查 probe/bait 文件里是否存在潜在重复。 |

## probe（19 个，含 Rust 版新增的 `easy-stampy`）

| 新命令 | 旧脚本名 | 说明 |
| --- | --- | --- |
| `probe easy-lastz` | `phyluce_probe_easy_lastz` | 用一条命令跑一次"简易" LASTZ 比对（一个文件对另一个文件）。 |
| `probe easy-stampy` | *（无，Rust 版新增）* | 用 [probebwa](https://github.com/GUIBA-EX/probebwa) 替代教程里手动调用的 `stampy.py`，一条命令跑完建索引、建哈希表、比对三步。 |
| `probe get-genome-sequences-from-bed` | `phyluce_probe_get_genome_sequences_from_bed` | 按 BED 文件里的坐标，从基因组中抽取对应 FASTA 序列。 |
| `probe get-locus-bed-from-lastz-files` | `phyluce_probe_get_locus_bed_from_lastz_files` | 把 bait 对基因组的 LASTZ 结果转成 **locus** 坐标的 BED 文件。 |
| `probe get-multi-fasta-table` | `phyluce_probe_get_multi_fasta_table` | 生成一张多路 FASTA 信息表。 |
| `probe get-multi-merge-table` | `phyluce_probe_get_multi_merge_table` | 生成一张多路合并信息表。 |
| `probe get-probe-bed-from-lastz-files` | `phyluce_probe_get_probe_bed_from_lastz_files` | 把 bait 对基因组的 LASTZ 结果转成 **bait** 坐标的 BED 文件。 |
| `probe get-screened-loci-by-proximity` | `phyluce_probe_get_screened_loci_by_proximity` | 对彼此距离过近的 bait，每组只保留 1 个 locus。 |
| `probe get-subsets-of-tiled-probes` | `phyluce_probe_get_subsets_of_tiled_probes` | 从多物种设计的 bait 集里，按需要的物种子集裁剪出对应 bait。 |
| `probe get-tiled-probe-from-multiple-inputs` | `phyluce_probe_get_tiled_probe_from_multiple_inputs` | 用多个输入基因组设计 bait。 |
| `probe get-tiled-probes` | `phyluce_probe_get_tiled_probes` | 用单个输入基因组设计 bait。 |
| `probe query-multi-fasta-table` | `phyluce_probe_query_multi_fasta_table` | 查询多路 FASTA 信息表。 |
| `probe query-multi-merge-table` | `phyluce_probe_query_multi_merge_table` | 查询多路合并信息表。 |
| `probe reconstruct-uce-from-probe` | `phyluce_probe_reconstruct_uce_from_probe` | 从 UCE bait 集反推设计时用到的 UCE locus 序列。 |
| `probe remove-duplicate-hits-from-probes-using-lastz` | `phyluce_probe_remove_duplicate_hits_from_probes_using_lastz` | 用 bait 自比自的 LASTZ 结果，筛掉跟其他 bait 匹配的疑似重复 probe。 |
| `probe remove-overlapping-probes-given-config` | `phyluce_probe_remove_overlapping_probes_given_config` | 按配置文件里的名单，从 bait 集里过滤掉指定 bait。 |
| `probe run-multiple-lastzs-sqlite` | `phyluce_probe_run_multiple_lastzs_sqlite` | 对多个基因组批量跑 LASTZ 搜索，结果写入 SQLite。 |
| `probe slice-sequence-from-genomes` | `phyluce_probe_slice_sequence_from_genomes` | 根据 `run-multiple-lastzs-sqlite` 的结果，从对应基因组切出匹配区域的序列。 |
| `probe strip-masked-loci-from-set` | `phyluce_probe_strip_masked_loci_from_set` | 从候选 bait 集里剔除掩蔽（masking）比例过高的序列。 |

## genetrees（5 个）

| 新命令 | 旧脚本名 | 说明 |
| --- | --- | --- |
| `genetrees generate-multilocus-bootstrap-count` | `phyluce_genetrees_generate_multilocus_bootstrap_count` | 配合按位点/按 locus 的 bootstrap 重抽样使用；现在不太推荐用这套流程。 |
| `genetrees get-mean-bootrep-support` | `phyluce_genetrees_get_mean_bootrep_support` | 给一批 gene tree，计算平均 bootstrap 支持率。 |
| `genetrees get-tree-counts` | `phyluce_genetrees_get_tree_counts` | 用对称差异比较一批 alignment 对应的 gene tree 拓扑，统计相同/不同拓扑的数量。 |
| `genetrees rename-tree-leaves` | `phyluce_genetrees_rename_tree_leaves` | 按配置文件里的新旧名字映射，重命名树的叶子节点；`--reroot` 可同时重新生根。 |
| `genetrees sort-multilocus-bootstraps` | `phyluce_genetrees_sort_multilocus_bootstraps` | 配合按位点/按 locus 的 bootstrap 重抽样使用；现在不太推荐用这套流程。 |

## ncbi（2 个）

| 新命令 | 旧脚本名 | 说明 |
| --- | --- | --- |
| `ncbi chunk-fasta-for-ncbi` | `phyluce_ncbi_chunk_fasta_for_ncbi` | 把一个 FASTA 文件拆成每份不超过 10000 条序列的分块（Sequin 文件的限制）。 |
| `ncbi prep-uce-align-files-for-ncbi` | `phyluce_ncbi_prep_uce_align_files_for_ncbi` | 把一批 alignment 整理成适合 `tbl2asn` 处理、用于提交 NCBI 的格式。 |

## utilities（8 个）

| 新命令 | 旧脚本名 | 说明 |
| --- | --- | --- |
| `utilities combine-reads` | `phyluce_utilities_combine_reads` | 按配置文件把同一批次的 reads 合并到一起。 |
| `utilities filter-bed-by-fasta` | `phyluce_utilities_filter_bed_by_fasta` | 用一份 UCE FASTA 文件过滤对应的 BED 文件。 |
| `utilities get-bed-from-fasta` | `phyluce_utilities_get_bed_from_fasta` | 给一个 bait FASTA 文件，生成对应位置的 BED 文件。 |
| `utilities merge-multiple-gzip-files` | `phyluce_utilities_merge_multiple_gzip_files` | 把同一个样本的多个 gzip 文件合并成一个；`--trimmed` 时按 R1/R2/singleton 分别合并已修剪的 reads。 |
| `utilities merge-next-seq-gzip-files` | `phyluce_utilities_merge_next_seq_gzip_files` | 合并 NextSeq 测序仪产出的多个 fastq gzip 文件（同一样本有时会拆成 4 个文件）。 |
| `utilities replace-many-links` | `phyluce_utilities_replace_many_links` | 用配置文件批量重写一批 symlink。 |
| `utilities sample-reads-from-files` | `phyluce_utilities_sample_reads_from_files` | 按比例从 fastq 文件中随机抽样 reads；Rust 版依赖 `seqkit`（原版依赖 `seqtk`），细节见 [已知差异](../rust-command-compatibility.md)。 |
| `utilities unmix-fasta-reads` | `phyluce_utilities_unmix_fasta_reads` | 把交错（interleaved）的 fastq 文件拆分成 R1、R2 和 singleton 三份。 |

## workflow（1 个）

| 新命令 | 旧脚本名 | 说明 |
| --- | --- | --- |
| `workflow` | `phyluce_workflow` | 运行各类 Snakemake 工作流的统一入口。 |
