# phyluce Rust 文档

`phyluce` Rust CLI 为 UCE 及一般系统发育基因组学流程提供统一命令入口。它覆盖
组装封装、contig 与探针匹配、UCE 序列提取、按位点比对、修剪、过滤、拼接、
统计，以及 probe、gene-tree 和常用数据整理工具。

本目录只说明 Rust 实现。原版 PHYLUCE 文档用于理解 UCE 工作流和分析背景；
安装方法、命令名称和兼容边界应以本项目文档及 `phyluce --help` 为准。

## 从这里开始

| 读者 | 推荐阅读顺序 |
| --- | --- |
| 第一次运行 | [构建与安装](user-manual-zh.md#2-构建与安装) -> [配置外部程序](user-manual-zh.md#3-配置外部程序) -> [端到端示例](user-manual-zh.md#17-端到端最小示例) |
| 从原版脚本迁移 | [命令形式](user-manual-zh.md#1-命令形式) -> [旧脚本名兼容](user-manual-zh.md#5-原版脚本名兼容) -> [已知差异](user-manual-zh.md#14-已知差异) |
| 查询具体命令 | [UCE 主流程](user-manual-zh.md#6-uce-主流程) -> [Alignment 工作流](user-manual-zh.md#7-alignment-工作流) -> 其他命令域 |
| 开发和验证 | [项目 README](../README.md#开发检查) -> [兼容性测试](user-manual-zh.md#16-兼容性测试) |

## 文档地图

- [README](../README.md)：项目结构、主要实现变化和开发检查。
- [中文用户手册](user-manual-zh.md)：安装、配置、完整流程、全部命令域、
  故障排查和复现建议。
- [旧命令兼容说明](../rust-command-compatibility.md)：旧脚本名映射、日志和
  迁移边界。
- [许可证](../LICENSE)：本项目的软件许可条款。

## 使用边界

- Rust CLI 不负责原始 reads 的接头去除和低质量碱基清理；进入组装前应使用
  合适的质控工具完成这些步骤。
- MAFFT、LASTZ、SPAdes、Gblocks、trimAl、RAxML、Snakemake 等仍是外部
  程序，只有调用相应命令时才需要安装并配置。
- Rust CLI 可准备拼接矩阵、分区文件和 gene-tree 辅助数据，但不替代最终的
  系统发育推断软件。
- 正式分析前，应先用小型数据验证命名规则、外部工具版本和输出目录结构。

## 上游参考

以下原版 PHYLUCE 页面适合了解方法背景和经典 UCE 工作流：

- [PHYLUCE 官方文档首页](https://phyluce.readthedocs.io/en/latest/index.html)
- [项目用途](https://phyluce.readthedocs.io/en/latest/purpose.html)
- [Tutorial I：UCE phylogenomics](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-1.html)
- [原版安装说明](https://phyluce.readthedocs.io/en/latest/installation.html)
- [引用说明](https://phyluce.readthedocs.io/en/latest/citing.html)

这些链接描述原版 Python 项目。Rust 实现的可用参数和输出行为以当前版本帮助
信息、兼容性测试和本地手册为准。
