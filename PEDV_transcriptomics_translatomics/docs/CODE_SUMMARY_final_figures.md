# 终稿图表（Figure 1–6、Sup Figure 1–19）数据分析代码总结

> 范围：`<PEDV_project>\final figure\final\新建文件夹\` 下 25 张终稿图（Figure 1–6，Sup Figure 1–19）。
> 查找范围：`<PEDV_project>\`（主要在 `batch3_agent_tasks\`）与 `<CoV_frameshift_project>\`，共扫描 1,406 个脚本。
> 路径缩写：**P** = `<PEDV_project>\batch3_agent_tasks\`；**K** = `<CoV_frameshift_project>\`；**FR** = `<PEDV_project>\final_figure_reproducibility\`（终稿数据包，是脚本的拷贝，不是原始位置）。
> 配套索引：`FIGURE_SCRIPT_INDEX.tsv`（逐面板：图 → 脚本 → 数据 → 状态）。

## 0. 怎么读这份总结

**可信度标记**（每张图都标了）：

| 标记 | 含义 |
|---|---|
| ✓ | 终稿图上的文字/数值、脚本里的字符串与 docstring、源数据三者对上了 |
| ◐ | 脚本与图内容高度吻合，但没有直接字符串证据，或终稿与脚本输出存在差异，需要你核对 |
| ✗ | 终稿图不是由本项目代码产生（湿实验、Prism、示意图、外部程序） |

**重要前提**：`.ai` 文件的置入路径全部被压缩，只有 Sup 16 保留了 `Supplementary_Figure_L_REVISED_v3_p2.pdf`（直接证据）。其余映射靠"读图内容 + 在脚本里搜图上的标签字符串 + 对照 DELIVERY/METHODS 文档"得出。**图上若有 Illustrator 手工改动（改名、删面板、改标注），代码输出与终稿会不一致，见第 6 节。**

**全项目统一口径**（下面各图不再重复）：
- 参考：`refdata_GIb_v1.1`，PEDV_GIb 全长 **27,959 nt**，1-based；宿主 Vero 为 ChlSab1.1（PK1 数据用猪 Sscrofa11.1）。
- 毒株统一为 GI-b（≈KP162057.1）；RNA-seq 为 **2 MOI（0.1/1）× 2 处理（CHX/HAR）× 5 hpi（3/5/8/12/16）每格 n=1**，Ribo-seq 一律 1 MOI；**零生物学重复**，各图不画误差棒、不做时间点显著性检验。
- 经典 ORF 坐标：ORF1a 293–12622 / ORF1b 12592–20613 / S 20610–24758 / ORF3 24758–25033（截短 91 aa）/ E 25361–25600 / M 25608–26288 / N 26300–27625。

---

## 1. 总览：图 → 主脚本

| 终稿图 | 内容 | 主绘图/分析脚本 | 标记 |
|---|---|---|---|
| **Fig 1** | A 设计；B reads 组成；C 全基因组 RNA 覆盖；D 病毒 RNA 动力学；E 3′/5′ 比；F 区域积累 | P`Fig1_MOI_regrouped\scripts\plot_Fig1_MOI.py`（B/C/D/F）、`plot_Fig1E_MOI.py`（E）；数据层 `fig1_moi_data.py` | A ✗；B–F ✓ |
| **Fig 2** | junction 景观 A–G | P`GIb_rnaseq_reanalysis\scripts\fig2_classify.py` → `fig2_plot.py` | ✓ |
| **Fig 3** | TRS-L 相似度 / 微同源 / 配对 A–H | P`GIb_rnaseq_reanalysis\figures\Fig3_work\`（`s3_*`、`s4_*`、`s5_*`、`metrics_*.py`） | ◐ |
| **Fig 4** | Ribo-seq QC A–E | P`fig4_v2_host_criteria\scripts\make_fig4_v2.py` | ✓ |
| **Fig 5** | gRNA 5′/3′ 翻译景观 A–D | P`Fig5_recolor_Fig4palette\plot_Fig5_recolored.py`；分析 P`fig5_terminal_ORF\scripts\`、P`fig5_gRNA5prime_figure_package\scripts\` | ✓ |
| **Fig 6** | A–D 序列扫描；E–F 起始/ORF 证据；G 质谱 | A–C：K`Result6_genomics_figures\scripts\figs_abc.py`；D：`figs.py`；E–F：P`HARTIS_v3\G4_Figure6_integration\scripts\`；G：`fig6.pzfx` | A–F ◐，G ✗ |
| Sup 1 | Western + RT-qPCR | `Sup Figure 1.pzfx`（GraphPad Prism） | ✗ |
| Sup 2 | RNA-seq 比对统计 / 相关热图 / PCA / reads 组成 | P`Fig1_MOI_regrouped\scripts\plot_SuppFig1_QC_MOI.py`（调用 `GIb_rnaseq_reanalysis\scripts\plot_SuppFig1C/D/E.py`） | ✓ |
| Sup 3 | 非经典 junction 的 5′/3′ 断点分布 | P`Fig2_Supp_junction_dynamics\scripts\s9_panelE_standalone.py` | ✓ |
| Sup 4 | 四类 junction：Vero / 猪体内 / LLC-PK1 | P`public_data_replication\ESV\biological_validation_v1\scripts\s06_figs.py`（A、B）、`s17_nc_corr_figs.py`（C、D） | ◐ |
| Sup 5 | 经典 sgRNA 精确断点变体 | FR`Sup_\plot_canonical_variants.py`（外部生成，见 5.5） | ◐ |
| Sup 6 | junction 类别组成 / 负担 / 归一化丰度（按 MOI） | P`Fig2_Supp_junction_dynamics\scripts\s4_plot_supp_figure.py`（A）、`s8_panelC_standalone.py`（B、C） | ✓ |
| Sup 7 | 30 个非经典 junction 簇行 z-score 热图 | 同目录 `s4_plot_supp_figure.py` + `_patch_layout*.py` | ✓ |
| Sup 8 | 非经典 sgRNA 推定 ORF 景观 + Top-20 | P`Noncanonical_ORF_classes\figures\make_noncanonical_figures_v2.py`、`make_top20_single_panel.py` | ◐ |
| Sup 9 | Top-20 junction 跨 Vero / 体内 / LLC-PK1 | P`public_data_replication\ESV\biological_validation_v1\scripts\s03_highfreq.py` + `s06_figs.py` | ◐ |
| Sup 10 | Ribo reads 去向、链向、反义、长度 | P`supfig12_v2_ABC\scripts\make_supfig12_v2_ABC.py`（A、D、E）；P`Supplementary_Figure_13\code\plot_supplementary_figure_13.py`（B、C） | ◐ |
| Sup 11 | 宿主 CDS 读框 × 读长（6 个库） | P`fig4_v2_host_criteria\scripts\make_grouped_frame_bars_v2.py` ← `read_length_frame_periodicity\scripts\01–05` | ✓ |
| Sup 12 | 病毒读框周期性 | P`Supplementary_Figure_13\code\plot_supplementary_figure_13.py` ← `prepare_panels.py` | ✓ |
| Sup 13 | gRNA 5′/3′ 候选位点按读框着色 | P`Fig5_3UTR_translation_audit\scripts\landscape_5p_style.py`、`landscape_3p_style.py` | ◐ |
| Sup 14 | 多序列比对的三读框 per-sequence 热图 | 用户自己的 FRASynPlot 流程：`<PEDV_project>\PEDV_ribo\`（`synplot2_vis.R`） | ✗（外部/早期） |
| Sup 15 | 翻译效率（Ribo/RNA）vs 长 ORF 容量 | K`Result6_genomics_figures\scripts\fig_E_onesided.py` | ✓ |
| Sup 16 | 起始密码子状态与家族 ORF 架构（按谱系） | K`Supplementary_Figure_L\scripts\f3_panelC.py`、`f4_make_figure_v3.py` | ✓（.ai 直接证据） |
| Sup 17 | N-x1 / N-x2 / ORF1b-x13 / S-x1 局部 P-site | P`Supplementary_Figure_N\PRICE29_locus\code\price29_locus_2x2_by_sample.py` | ✓ |
| Sup 18 | LLC-PK1 与 Vero 的 P-site 对照 | P`PK1_VERO_ORF_VALIDATION\scripts\corr_figs_editable.py`（A）、`nx_local_16h_v11_editable.py`（B） | ◐ |
| Sup 19 | N-x1 PRM 色谱 + 肽段位置 | P`PRM_Nx1_redraw\scripts\05_redraw.py`（曲线为数字化近似） | ◐ |

---

## 2. 公共流程模块（图里引用的"上游分析"）

### M1　RNA-seq：比对 → junction → 分类（Fig 1–3、Sup 2–9）
| 步骤 | 脚本 | 关键点 |
|---|---|---|
| 比对 | P`GIb_rnaseq_reanalysis\scripts\run_align_GIb.sh` | STAR 2.7.11b（`merged\star_ChlSab1.1_PEDV_GIb_oh149`），junction 友好参数：`--alignIntronMax 30000 --alignMatesGapMax 30000 --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 --outFilterMismatchNmax 6 --outFilterMismatchNoverLmax 0.08 --outFilterMultimapNmax 10 --alignSJoverhangMin 8`；`samtools markdup`（坐标法去重，仅供覆盖度）；`samtools depth -a` 出逐碱基覆盖 |
| 独立 junction 比对 | 同上 | **minimap2 `-ax splice -un -G 30000 -k 13 -w 5 --secondary=no`**，病毒单独索引，未去重原始 reads；**所有 junction 终稿数字来自 minimap2，不是 STAR**（STAR SJ.out.tab 一致性未检验，见 Fig 3 方法文档的限制项） |
| junction 提取 | `extract_junctions_GIb.py` | 逐 gap 独立计算（修复多 gap read 共享 effective_leader_start 的 bug）；`MIN_GAP=50`；输出全部原始候选 |
| 分类 | `fig2_classify.py`（Sup 3–9 的 `s1_core_counts.py` 原样沿用） | 保留 primary、非重复、`NH==1`；`MAPQ≥20` 且两侧锚长 **≥15**（leader 锚 = gap 前所有 M 之和）；同一 fragment 跨同一 (donor, acceptor) 只计 1 次；**TRS-L 依赖 = donor 在 leader 区 [55,72]**（核心 CTAAAC 67–72）；**canonical = acceptor 在 S 20597 / ORF3 24718 / E 25342 / M 25597 / N 26285 ±5**；TRS-L 非依赖按 gap=acc−donor−1 分 **local <1 kb / distant ≥1 kb**；非经典 junction 的受体 ORF 后果由 `first_atg()` + `cds_of()` 判 in-frame / frameshift / 无 AUG / 非编码 |
| 产物 | — | 43,753 个事件，2,674,258 条 junction-spanning fragments；canonical 95.55%、TRS-L-dep nc 1.22%、distant 1.72%、local 1.52%（Fig 2B） |
| atlas | `build_sgRNA_atlas_GIb.py` | 23 样本合并链式聚类（相邻 ≤5 nt）、peak-aware 拆分（簇宽 >12 nt 时按读数谷值 ≤较低峰 50% 拆分）、置信分级 |

### M2　Ribo-seq：预处理 → 比对 → P-site → 读长集合（Fig 4、Sup 10–13、Fig 5/6 的底层数据）
| 步骤 | 位置 | 关键点 |
|---|---|---|
| 预处理 | 服务器 `analysis_batch23/Riboseq_qc_mapping/scripts/run_one.sh`（**本地没有该脚本，只有 P`batch23_riboseq_qc_mapping_output\make_report.py` 汇总**） | cutadapt（接头 RA3，`--discard-untrimmed`）→ bowtie2 rRNA/ncRNA 扣除 → STAR；核对 bowtie2 输入 = cutadapt 输出、STAR 输入 = bowtie2 未比对 reads |
| 比对 | P`riboseq_GIb_mapping_psite\scripts\server\run_star_one.sh` | STAR 2.7.11b，`refdata_GIb_v1.1\merged\star_ChlSab1.1_PEDV_GIb_oh40`；**`--outFilterMultimapNmax 1 --outFilterMismatchNmax 2`**（所以"多位点"全部落在 unmapped/too-many-loci，不会出现 NH>1） |
| BAM 统计 | `gib_bam_stats.py` | PEDV reads **一条不删**；读长 = 比对区块长度（CIGAR M/=/X/D 之和，跳过 N，不含 soft-clip）；输出逐读长计数、5′ 端矩阵、重复率 |
| 宿主 P-site | `host_psite_pass.py` | Ensembl canonical protein_coding；起始/终止 ±60 nt metagene、逐基因 CDS 三读框计数；P-site = 5′ 端 + offset |
| offset 与读长集合 | `build_B_offset_sets.py`（Module 1 方法 B，`module1_rpf_qc\scripts\m1_build_tables.py` 为母版） | **只用宿主起始密码子校准 offset**（PEDV 起始位点因 leader soft-clip 不参与）；C1 稳定 P-site（≥500 reads、峰占比 ≥10%、峰 ≥次高 1.2×、两半稳定）、C2 3-nt 周期性（ORF 唯一区两端各去 45 nt，≥1,000 条、最高帧 ≥次高帧 1.5×）、C3 frame0 ≥0.40；四套读长集合 `strict_core`（29/30）/ `fixed_29_30` / `expanded` / `short_RPF`；≥35 nt 为 LPF（保留，不进密码子分辨率分析） |
| PEDV P-site 长表 | `build_C_psite_tables.py` | `psite_counts_PEDV_GIb.tsv.gz`（744,452 行：sample / strand / position / read_length / offset / length_category / genome_frame / orf_frame） |
| 共享计算层 | P`shared_riboseq_metrics\scripts\build_shared_metrics.py` | Sup 10–12 共用；每个分母恒等式都断言（non_rRNA == STAR input；host+PEDV == uniquely mapped；七类之和 == raw） |

**口径说明（影响多张图）：** CHX 16 h 在 Fig 4 / Sup 10–13 / Sup 17 等已用 **XN26040813**；早期的 Fig 4 v1（`fig4_riboseq_qc\scripts\make_fig4.py`）和扰动筛选等仍有 065，**当前终稿读长集合是各库"通过宿主判据的读长"（host_criteria：offset high、≥5,000 CDS-body P-site、F0≥0.5、最高帧/次高帧≥1.5），不是 strict_core**；Fig 6 E–F 的 PRICE/HARTIS 是 **29 nt** 口径。

### M3　ORF calling（Fig 6E–F、Sup 17 的候选来源）
- PRICE 1.0.5（Gedi 1.0.6f 源码编译）、Ribo-TISH 0.2.8、ribotricer 1.5.0、RiboCode；脚本在 P`orf_calling_GIb\scripts\`（及 `scripts\server\`、`price_813\`、`ribotish_29nt\`、`integrate\`）与 P`orf_calling_GIb\_price_final\_29nt_work\scripts\`。
- **Fig 6 用的是 PRICE 29 nt**：宿主训练**不加 `-filter`**，病毒 `-filter 29:29`，宿主模型复用到 + / − 链（`p29all_run.sh`；4 个范围：060、811、813、merged_CHX）。
- 项目内的已知坑：PRICE 在病毒上起始预测不可靠（常给出 ATT）；宿主模型不稳定；Ribo-TISH 的默认 BH 是各运行内部校正，跨运行 q 不可比。

### M4　HARTIS_v3（Fig 6E–F 的 HAR 起始证据）
- 位置 P`HARTIS_v3\`（`scripts\v3_01–v3_06`、`G2b_host_background\`、`G3_pedv_initiation\`、`G4_Figure6_integration\`）。
- 思路：HAR 起始近端窗（W1）的 +1 框占据，对照**宿主位点级条件经验零分布**（5′ 端序列上下文、距起止作协变量）；`in_support`（EMP_Zs）与 `capped_reference`（探索性，未经病毒高深度校准）**分开画，不合并**；同一文库多个读长只计一次文库；"最小名义 P"为未校准、非联合 P。
- 预注册文档：`PREREGISTRATION_HARTIS_v3*.md`、`G2_FINAL_HARTIS_v3.md`、`BLINDING_LOG.md`。

### M5　比较基因组学（Fig 6A–D、Sup 14–16）
- 序列库：1,083 条 PEDV 基因组（来源 K`N1_conservation_analysis\data_external\` 与 `results\PEDV_genomes_all.fasta`，经 Sup L 审计）；同源坐标映射 K`Supplementary_Figure_L\scripts\mapper.py`、`p1_map_all.py`；谱系分型 `p2_typing_dedup.py`（S 基因投影）。
- stop-anchored 候选池：K`Supplementary_Figure_L\scripts\p3a_reference_pool.py`（三读框、8 种起始密码子 = AUG + CUG/GUG/UUG/AUA/AUU/AUC/ACG、≥8 aa、完整 ORF；归属规则 fixed assignment）；保留性 `p3b_retention.py`；独立重扫 `p4a_independent_scan.py`；同义替换零模型 `p5a_synonymous_null.py`。
- Result 6 的统计 K`Result6_genomics_figures\scripts\`：`g1_null.py`（组成零模型：保持氨基酸序列、按各 ORF 自身密码子使用重抽同义密码子，10,000 次，seed 20261019）、`g1_conservation.py`、`g2_length.py`、`g3_g4_tables.py`、`g4_activity.py`、`g5_null.py`/`g5_null_het.py`/`g5_test.py`（Poisson-binomial 精确尾 + BH）、`g5_synplot_user2.py`（SynPlot2 式检验，用户的去冗余基因比对）。

---

## 3. Figure 1–6 逐图说明

### Figure 1　RNA-seq 概览
- **分析问题**：感染进程中病毒 reads 占比、基因组覆盖形态、总病毒 RNA 动力学、3′/5′ 覆盖比、区域积累。
- **脚本**：数据层 P`Fig1_MOI_regrouped\scripts\fig1_moi_data.py`（复用 `plot_Fig1DEF.py`、`plot_RNA00_read_composition.py`、`plot_figures_GIb.py` 的口径）；B/C/D/F：`plot_Fig1_MOI.py`；E：`plot_Fig1E_MOI.py`。
- **口径**：
  - B：PEDV / Host / Other 占比，分母 = 全部 QC 后 read pairs（`count_read_composition.py`：以 STAR 主/次比对涉及的 contig 类别归类为 pedv_only / host_only / cross / unmapped）。终稿 B 的 16 hpi PEDV 约 31%，与 Sup 2E 的 MOI=1 一致，**推测 B 取的是 MOI=1 子集（`*_main_MOI1*`）**。
  - C：逐碱基深度 / 总深度 ×10⁶ 的 101-nt 滑动平均，按 MOI 分上下两格；CHX 实线、HAR 虚线。
  - D：唯一比对（MAPQ=255）PEDV fragment 占映射 reads 的 %，**未去重（raw）为主**。
  - E：3′/5′ 覆盖比，窗口 `3p_main`/`5p_main`，raw，MAPQ=255。
  - F：区域（Leader、ORF1a、ORF1b、S、ORF3/E/M、N/3′UTR）相对覆盖热图，色阶 ±4，格内不标数字。
  - **无误差棒、不跨 MOI 平均**（原图把两个 MOI 当 n=2 画 SEM，是已更正的错误）。
- **Panel A**：示意图，无代码。图中仍写 "MOI = 0.1 or 1"，与更正后的设计一致。

### Figure 2　junction 景观
- **脚本**：`fig2_classify.py`（分类，见 M1）→ `fig2_plot.py`（`panel_A` … `panel_DG`）；数据包在 FR`Figure_2\`。
- **A**：5′ × 3′ 断点二维直方（100-nt bin，log₂，截断于 2¹⁶）+ 边缘直方；橙色标记 a–f 为热点（`hot_list()`，源数据 `Fig2A_hotspots_source_data.csv`）。
- **B**：四类 junction 的 reads 占比（n = 2,674,258）。**C**：canonical 各 sgRNA 的 reads（对数）+ TRS-L 依赖非经典物种（n=376）按 reads 排序，颜色 = 受体 ORF 后果。
- **D–G**：四类 junction 的弧线图；弧线颜色按 frame（in-frame / out / 无 AUG / 非编码）；标题里的比例（如 "in-frame 0.10 / out-of-frame 0.86"）来自 `class_pct()`/`stats()`。
- **敏感性（已在 `Fig2_DELIVERY.md` 记录）**：leader 区 57–67 → 1.08%、50–75 → 1.26%；distant 阈值 500 nt 1.99% / 1 kb 1.72% / 2 kb 1.35% / 10 kb 0.96%（`events_*.tsv`）。

### Figure 3　TRS-L 相似度、微同源与配对
- **脚本**：P`GIb_rnaseq_reanalysis\figures\Fig3_work\`：`s3_lib.py`（`mh_vec`、`chain`、`boot_mean`、`p2`、`holm`）、`metrics_sim_mh.py`（微同源、二核苷酸打乱 `dinuc_shuffle`、`best_core`）、`metrics_mfe.py`（RNAduplex 能量窗口）、`s3_panels.py`/`s4_panels.py`/`s5_panels.py`（各轮绘图）、`step1_join_check.py`（事件表 69,208 → 43,753 的对账）。
- **统计方法（`Fig3_statistics_methods.md`，逐条核对过）**：
  - 全部随机化 seed 42，重复 10,000 次，经验双侧 P = (1+#)/(1+N)，最小可得 P = 9.999e-5；**分析单位是唯一 (donor, acceptor) junction，reads 不作权重**。
  - 断点约定：1-based、**左对齐**（donor = 最后一个保留的 leader 碱基，acceptor = 第一个 body 碱基）；左滑 L=0 的事件 43,752/43,753，唯一例外 59→25344（canonical，L=2）；故"歧义窗口 = 微同源 + 1"是构造恒等式。断点只来自 minimap2。
  - TRS-L 相似度：41 个窗口偏移（−20…+20）取最大匹配数（背景用同样的最大化，背景均值 6.6/12 而非 3.2）；背景 = 二核苷酸打乱的 52-nt 受体上下文，每个受体 200 次。
  - 微同源：背景 = 类内 pair 置换（swap chain，2 sweeps/记录样本，burn-in 20），distant 保持 gap ≥1,000 nt，local 保持 0<gap<1,000 nt；canonical 因只有 5 个受体、donor 近乎固定，改用位置背景。
  - MFE：ViennaRNA RNAduplex，donor/acceptor 各 30 nt 窗口；(a) donor 窗口反向互补对受体窗口，(b) 两条正链窗口直接配对；RNAcofold 交叉核对 ρ=0.88（只证算法稳健，不证能量模型）。
  - 多重检验：Holm，按 family（4 个组或 3 个类）。
- **◐ 需要你核对**：
  1. 终稿 3D 标注的微同源增量是 **+5.1 nt (7.8×) / +1.07 nt (2.2×) / +0.75 nt (1.9×) / +0.96 nt (2.0×)**，而 `Fig3_statistics_methods.md` 里写的是 +5.05 / +1.09 / +0.771 / +1.06，**Dep/Dist/Loc 三个数字略有差异**（local 差 0.1 nt），说明终稿是更晚一次重跑的结果，或手工标注；
  2. `Fig3_work` 里 `s3_panels.py` 只有 A–G 七个 draw 函数，`s5_panels.py` 为 A–E，而终稿是 **A–H 八个面板**；**没有找到直接生成 8 面板拼图的单一脚本**，推测是不同轮次的面板在 Illustrator 里组合。

### Figure 4　Ribo-seq 文库 QC
- **脚本**：P`fig4_v2_host_criteria\scripts\make_fig4_v2.py`（`draw_A`…`draw_E`、`combined`、`single`）。**样本 = CHX 5/8 h（XN26080060、XN26040811）、CHX 16 h = XN26040813、HAR 5/8/16 h（XN26080061/62/64）**；读长集合 = host_criteria（25–35 nt）。数据全部取自 RGP_v1.3 与服务器一遍统计，不重算 offset（813 的 offset 由 Module 1 方法 B 判定函数在本地逐行复刻）。
- **脚本内部字母与终稿字母不同（已核对函数）**：脚本 `draw_D`（`exfrac`：宿主 P-site 在 CDS / 5′UTR / 3′UTR / other exon-ncRNA 中的占比）→ **终稿 A**；脚本 `draw_A`（读长分布，x 轴 18–50 nt，"Fraction of reads"）→ **终稿 B**；脚本 `draw_B`（宿主 CDS-body 三读框占比）→ **终稿 C**（红点 = 单个文库，RNA-seq 伪 P-site 作阴性对照，`rna_pseudo_psite.py`：5′ 端 +12 nt）；脚本 `draw_C`（起始/终止密码子 metagene，每 10⁶ 个 CDS P-site）→ **终稿 D**；脚本 `draw_E`（6 个库 PEDV 全基因组 P-site 轨道，对数）→ **终稿 E**。即终稿在 Illustrator 里按新顺序重排过；对账时按函数名，不要按字母。
- **终稿 E 的最大值**（CHX 5/8/16 h = 6.7 / 10.6 / 11.7 ×10⁴，HAR 5 h = 2.6）对应 **"全部读长"版本**（README 说 HAR 5 h 选定读长后最大值降为 1.0×10⁴），即 `Fig4_v2_combined_E_all_lengths` 的口径。
- CDS-body 定义：CDS 两端各去 45 nt；终止窗口读框 = (rel + CDS 长度) mod 3。

### Figure 5　gRNA 5′ / 3′ 翻译景观
- **绘图**：P`Fig5_recolor_Fig4palette\plot_Fig5_recolored.py`（`--base` 指向 `fig5_gRNA5prime_figure_package`；由 GPT 提供的 `plot_Fig5_stylematched_*` 改配色而来，见 5.4）。
- **分析**：P`fig5_terminal_ORF\scripts\`（`step0_junction_remeasure.py` → `step1_build_catalog.py`（gRNA、S、ORF3、E、M、N 转录本模型 + 末端 ORF 枚举，仅序列）→ `step2_evidence.py`（HAR 起始窗 / CHX body / spike / RNA 背景，固定阈值描述符，**无综合评分、无二元调用**）→ `step3_spike_frame.py` → `step4_tables.py` → `step5_figures.py`）与 P`fig5_gRNA5prime_figure_package\scripts\g1_gRNA5p_analysis.py`、`g2_gRNA5p_tables.py`、`q1/q2/q3`（QC 与草图）。
- **口径**：RGP_v1.3 strict_core P-site（+ 链），RNA 深度对数；5′ 区把近同源起始（CUG、GUG、UUG、ACG、AUA、AUC、AUU）与 AUG 分开标注；"有信号"的候选红线、"无信号"蓝线。
- 3′ 区（C、D）：N 末 300 nt + 3′UTR 的候选 ORF 景观；审计脚本 P`Fig5_3UTR_translation_audit\scripts\audit3utr*.py`（`ANALYSIS_PLAN_predeclared.md` 预声明）。

### Figure 6　PEDV 基因组中的移码 ORF 筛选
- **A–D（序列层）**：K`Result6_genomics_figures\scripts\`：`stats_abc.py`（卡方检验）→ `figs_abc.py`（A、B、C）、`figs.py`（D）；数据 `g1_null.py` 的组成零模型。
  - A：每 kb CDS 的终止密码子 / AUG 在 0、+1、+2 框的数目（误差线为置信区间）；B：stop-anchored family 最大肽长直方图（5-aa bin），>40 aa：+1 框 **68**、+2 框 **3**；C：每 kb 长（>40 aa）移码 family 数（**+1 与 +2 合并**：ORF1a 32/12.33 kb=2.6、ORF1b 18/8.02=2.2、S 12/4.15=2.9、ORF3 1/0.276=3.6、M 3/0.681=4.4、N 5/1.326=3.8，我按坐标复算与图上数字一致）；D：+1 框 family 在基因组上的位置，紫色 = 含 AUG、金色 = 仅近同源起始，红框 = featured family。
- **E、F（证据层）**：P`HARTIS_v3\G4_Figure6_integration\scripts\`：`fig6_build_tables.py`（整合 HARTIS 起始占据 × PRICE 29 nt ORF 证据）→ `fig6AB_genome_combined.py`（E 上下两格：HAR HARTIS W1 的 −log₁₀P 与 PRICE(CHX 29 nt) 的 −log₁₀P）、`fig6C_HAR_PRICE_joint_v2.py`（F 联合散点；v3_capped 版是探索性）。每个 family 一个点；无有效 P 的不绘制（不赋 P=1）；E-x1 在 PRICE 中被记为经典 M 读框，故不画。
- **G**：PRM/质谱肽段强度柱形图，**来自 GraphPad Prism 文件 `final figure\final\fig6.pzfx`**（表：Control / PEDV 12h / PEDV 16h；行：ORF1a、ORF1b、S、ORF3、E、M、N、N1），**不是本项目代码产物**。
- **◐ 需要你核对（终稿与代码输出有差异）**：
  1. `figs_abc.py` 生成的 **C 排除了 ORF3 与 E**（标题 "ORF3, E excluded"，5 根柱），终稿 C 含 ORF3（6 根柱）；
  2. 代码的 B 是 "families + start isoforms" 的 2×2 面板并带 `(max N aa)`，终稿 B 只有 families 两格，且文字 **"max :" 后面是空的（数字缺失）**；
  3. 终稿 C、A 没有标注卡方 p 值（代码里有）。
  → 终稿 Fig 6A–D 经过 Illustrator 手工改动，或来自 `figs.py` 的另一版本；**以终稿为准重做源数据对账**。

---

## 4. Supplementary Figure 逐图说明

### Sup 1　湿实验（Western + RT-qPCR）
- A/B Western（N 与 β-actin，CHX/HAR × 5/8/12/16 h，Batch 1 = 0.1 MOI、Batch 2 = 1 MOI）为胶图；C–F RT-qPCR 柱形图来自 `Sup Figure 1.pzfx`（表：Batch1/2 × HAR/CHX × 3/5/8/12/16 h）。**无分析代码，原始 Ct/灰度值不在项目里。**

### Sup 2　RNA-seq QC
- P`Fig1_MOI_regrouped\scripts\plot_SuppFig1_QC_MOI.py`（复用 `GIb_rnaseq_reanalysis\scripts\plot_SuppFig1C.py`（A/B 文库与比对统计）、`plot_SuppFig1D.py`（C 样本相关，宿主 log₂(TPM+1) Pearson）、`plot_SuppFig1E.py`（D PCA，宿主基因 VST））。数据层 `compute_alignment_summary.py`、`count_read_composition.py`。
- 数值不重算，只把每个文库按 MOI × 处理 × hpi 逐一标注（旧图把同条件两个 MOI 当"两个重复"，并用灰线连接）；PCA 里 MOI 用空心/实心。E 面板为各条件 PEDV/Host/Other 占比。

### Sup 3　5′/3′ 断点分布（非经典 junction）
- P`Fig2_Supp_junction_dynamics\scripts\s9_panelE_standalone.py`。竖向 5 行 = 3/5/8/12/16 hpi，2 列 = MOI 0.1 | 1；x = 病毒基因组坐标，y = 每百万病毒 fragment 的 junction reads（对数）。**100-nt bin；同一 MOI × hpi 的 CHX + HAR reads 相加后除以两库病毒 fragment 之和；仅非经典（TRS-L 依赖 + 非依赖）；每条 junction read 在 5′ 与 3′ 直方图各贡献一次。**描述性，无检验。

### Sup 4　四类 junction：Vero / 猪体内 / LLC-PK1
- P`public_data_replication\ESV\biological_validation_v1\scripts\s06_figs.py`（`Fig1_four_class_junction_composition`：A、B），`s17_nc_corr_figs.py`（`FigA_scatter`：C、D）。分析：`gx_junction.py`/`s16_nc_correlation.py`（`gII_cross_strain_validation_v1\scripts\`）。
- **口径**：猪体内 = LJX 三只仔猪（P1/P2/P3，GII-b）、LLC-PK1 = 单样本（GII-b），**参考为各自共识参考（`LJX_merged_denovo` / `LLCPK_denovo`），acceptor 经 `asm20` 全基因组比对 lift 到 GI-b 坐标后再判类、再聚类（相邻 ≤5 nt）；不按原始坐标匹配不同毒株**；分类规则 `cls()` 原样来自 `p5_pe_junc.py`；深度不同，distinct junction 数不可直接比（报告里另给等深度）。
- C、D 的 Spearman（体内 0.24；LLC-PK 0.60）是跨宿主 nc 丰度相关。
- **◐**：图中标签（"In vivo P1…"）在 Illustrator 里改过，脚本里没有这些字符串，只能靠数字对账（`ESV_BIOLOGICAL_RESULTS.md` 表：canonical 95.55 / 96.54–97.43 / 85.16%）。

### Sup 5　经典 sgRNA 精确断点变体
- FR`Sup_\plot_canonical_variants.py`（输入 `Fig2D-G_arcs_source_data.csv`，输出 `canonical_exact_breakpoint_variants.csv`、`canonical_variant_summary.csv`）。A：每个 sgRNA 的 donor/acceptor 偏移散点（相对该类 top1 精确 junction，圆面积 ∝ reads）；B：top1/2/3/other 占比条（S 97.6%、ORF3 97.7%、E 97.2%、**M 42.4/31.5/25.8%**、N 99.6%）。**事件是精确 (donor, acceptor) 对，不是独特 RNA 分子。**这个脚本**只存在于数据包 `Sup_` 目录**，原始生成位置未找到。

### Sup 6　junction 类别组成与负担（按 MOI）
- A：P`Fig2_Supp_junction_dynamics\scripts\s4_plot_supp_figure.py`（上层分母 = 同库全部 classified junction reads，下层 = 仅非经典，**两层分母不同**）；B（Relative burden = 类 reads / 同库 classified junction reads）、C（Normalized abundance = 类 reads / (病毒 fragment ÷10⁶)）：`s8_panelC_standalone.py`。数据层 `s1_core_counts.py`（逐样本计数矩阵）→ `s3_panel_stats.py`。分母（病毒 fragment）= minimap2 primary 且 MAPQ≥20 的去重 read name，由服务器 `s2_server_denominator.sh` 统计。圆 = MOI 1，空心 = MOI 0.1；圆 = CHX、方 = HAR。

### Sup 7　30 个非经典 junction 簇的行 z-score 热图
- `s4_plot_supp_figure.py`（+ `_patch_layout*.py`）。行 = 非经典事件簇（≥100 pooled reads、在 ≥5 个库中各 ≥3 reads 检出）；值 = log₂(每百万病毒 fragment 的 reads + 1) 的**行内 z-score（色阶 ±2.5，只表示相对变化，不可比较事件间绝对丰度）**；斜线格 = 观察 0 且期望 <3（不可判定，不是真 0）。簇规则：TRS-L 依赖按 acceptor 链式聚类（≤5 nt），非依赖按 donor 与 acceptor 同时 ≤5 nt 单连接聚类；`DET_MIN=3`。†（27620→27682）注明 "Authenticity under audit"。

### Sup 8　非经典 sgRNA 推定 ORF 景观 + Top-20
- P`Noncanonical_reclass\reclassify.py`、`build_new_tables.py`（规则）→ P`Noncanonical_ORF_classes\figures\make_noncanonical_figures_v2.py`、`make_top20_single_panel.py`。**TRS-L 依赖非经典 junction（1,173 个）按其 sgRNA 能编码的 ORF 分类**：truncated（第一个下游 AUG 与受体所在 CDS 同框）/ frameshifted（异框）/ out-of-frame ORF（受体 CDS 内无 AUG 或非编码，但 sgRNA 含 ATG…stop ORF）；`translates_canonical`（sgRNA 含正常框的 ORF）**不计入**。ORF 身份 = 起始 ATG 位置，支持度 = 支持 junction 数。
- **命名冲突风险**：图里的 "N-x4、M-x2、M-x3、S-trunc、3UTR-1" 是**非经典 sgRNA 衍生 ORF 的名字**；Figure 6 里的 "N-x1…N-x5、ORF1a-x22" 是 **+1 框 stop-anchored family**。**两套名字互相重叠（两边都有 N-x4），指的不是同一个对象**，图注里需要区分。

### Sup 9　Top-20 junction 的跨宿主复现
- P`public_data_replication\ESV\biological_validation_v1\scripts\s03_highfreq.py`（逐样本 Top20、跨样本同源归并、与 Fig 2 簇对应——匹配区间两端各放宽 3 nt、TRS 相似度/微同源直接复用 `build_sgRNA_atlas_GIb.py` 的函数逻辑、右侧同源用 module5_9 engine 的 `junction_homology()`）+ `s06_figs.py`（`Fig4_top20_noncanonical_reads_reproducibility_conservation`）。A = Vero 20 样本 reads，B = 猪体内三只仔猪的 reads，C = LLC-PK1 reads；"20/20" 表示 20 个 Vero 库中检出的库数。◐（理由同 Sup 4）。

### Sup 10　Ribo-seq reads 去向、链向与读长特征
- A（分母链：raw → 接头/长度丢弃 → rRNA/ncRNA → 宿主唯一 / PEDV 唯一 / >1 位点 / unmapped；**注意 Nmax=1，所以"多位点"是 STAR 的 too-many-loci 类，不是真多比对为 0**）、D（中位读长）、E（29–31 nt 占比）：P`supfig12_v2_ABC\scripts\make_supfig12_v2_ABC.py`（取自 `shared_riboseq_metrics`）。
- B（PEDV 链向，含 strand-assignable reads 总数）、C（PEDV (−) reads % vs 宿主反义 P-site %）：P`Supplementary_Figure_13\code\plot_supplementary_figure_13.py` 的 panel A，数据 `prepare_panels.py`（`host_antisense_fraction_of_host_psites`）。
- 链向判读注意：Ribo-seq 文库保留链信息（RA3 接头连接），但 `prepare_panels.py` 的注释写明"未找到建库记录"，链特异性只是被支持而非被证实；负链占比 0.3–5.9%（同为 8 hpi 的 CHX 与 HAR 相差约 17 倍）。这里**不对负链做生物学解释**（我在 Sup 13 任务里预测它是技术产物，该预测是否已被检验我没有核对到结果）。

### Sup 11　宿主 CDS 读框 × 读长
- P`fig4_v2_host_criteria\scripts\make_grouped_frame_bars_v2.py`，数据 P`read_length_frame_periodicity\scripts\01_frame_counts.py`（宿主逐读长三读框，CDS-body）、`02_frame_ambiguity.py`、`03_tables_audit_figures.py`、`05_extra_sample_XN26040813.py`（813 同法复刻 offset/QC）。红点 = 通过宿主判据（offset high）的读长；虚线 = 1/3 无偏好基线。

### Sup 12　PEDV 读框周期性
- P`Supplementary_Figure_13\code\prepare_panels.py` → `plot_supplementary_figure_13.py`。ORF 唯一区 = [起点+45, 终点−45] 去掉与其他 ORF 重叠的碱基；frame = (P-site − ORF 起点) mod 3。主口径 host_criteria，另给 strict_core / 29–31 nt / all_accepted。C 热图的 **f = (F0 − 1/3)/(0.80 − 1/3)**，valid ≥1,000 reads、low depth 100–999（斜线）、missing <100。
- **注意（已在项目内确认）**：f 的上限 0.80 在 811/813/064 的 strict_core 上被突破（f>1），所以热图**不能当定量轴**，ORF1a/ORF1b 在 5/6 个库里 frame0 落在 1/3 本底（非核糖体信号）。

### Sup 13　gRNA 5′/3′ 候选位点按读框着色
- P`Fig5_3UTR_translation_audit\scripts\landscape_5p_style.py`（5′ 三窗：3–45、75–135、255–325；读框相对 ORF1a，青绿 = ORF1a 框、橙 = +1、紫 = +2）、`landscape_3p_style.py`（3′ 三窗；读框相对被标注候选自身，青绿 = 候选读框）；`candidate_profiles.py` 为逐候选小图。P-site = 正式 6 库（CHX 16 h = 813）。◐（面板命名对应基于 docstring）。

### Sup 14　多序列比对的三读框 per-sequence 热图
- **非本项目代码**。来自用户的早期 FRASynPlot（SynPlot2）流程：目录 `<PEDV_project>\PEDV_ribo\`；`流程.txt` 记录的步骤：头部对齐 → **MAFFT 多序列比对** → `cd-hit-est -c 0.99` 去冗余（`pedv_nr.fasta`）→ 提取 ATG ORF → `synplot2`；ORF 提取用 tBLASTn（e-value 1e-5）+ bedtools merge；热图由 `ORF1a\synplot2_vis.R`（ggplot2 + patchwork，`{prefix}_main.pdf` = synplot + per-sequence codon heatmap，3 个读框）画出。图中 7 个区域对应 `ORF1a/ORF1ab/spike/ORF3/envelope/membrane/nucleocapsid` 各自目录下的 `frasynplot_main.pdf`。
- 图上 "Mnvelope"（应为 Membrane）是拼写错误，需要在终稿修正。

### Sup 15　翻译效率 vs 长 ORF 容量
- K`Result6_genomics_figures\scripts\fig_E_onesided.py`：**TE = Ribo P-site 密度 / RNA 相对深度**（Ribo：`Fig6_canonical_background_audit.tsv` 的 CHX/HAR P_total，strict_core 29/30 nt、offset +12；RNA：`GIb_rnaseq_reanalysis\coverage\*.depth.tsv`，23 个库、非链特异），y = >40 aa +1/+2 family 每 kb 数；**n=6（ORF1a、ORF1b、S、ORF3、M、N；E 排除），单侧精确置换检验（720 次排列）**，只对"TE 越高、长 ORF 越多"这一方向假设。CHX ρ=0.83（p=0.029），HAR ρ=0.09（p=0.46）。**n=6 功效极低，且 N、M 的高覆盖含嵌套 sgRNA 的贡献，不等于翻译效率更高。**

### Sup 16　起始密码子状态与家族 ORF 架构（按谱系）
- K`Supplementary_Figure_L\scripts\f3_panelC.py`（C1：families × lineages 的堆叠条）、`f4_make_figure_v3.py`（v3 版式，`Supplementary_Figure_L_REVISED_v3_p2.pdf` = 此图，**.ai 内的直接证据**）。数据链：`p1_map_all.py`（同源映射）→ `p2_typing_dedup.py`（谱系分型）→ `p3a_reference_pool.py`（候选池）→ `p3b_retention.py`（保留性）→ `p5a_synonymous_null.py` → `p6_aggregate.py` → `p7_panels_and_qc.py` → `r2_audit_tables.py`。
- 谱系样本量 GI-a 23 / GI-b 30 / GII-a 752 / GII-b 96 / GII-c 158 / All 1,083。**这套分型来自 Sup L 的重新分型（S 基因投影），与早期 `genotype_normalized.tsv`（GIa=0、GIb=168）不同**，两套不能混用。起始状态六类：exact non-AUG 保留 / AUG 保留 / 近同源保留（不同密码子）/ 起始类别切换 / 丢失合格起始 / 不可评估（不可评估不进分母）。

### Sup 17　四个候选的局部 P-site
- P`Supplementary_Figure_N\PRICE29_locus\code\price29_locus_2x2_by_sample.py`（单候选版 `price29_locus_figs.py`、`price29_locus_1x4_by_sample.py`）。每个候选一格：ORF 示意 + 6 条 P-site 轨道（CHX 5/8/16 h、HAR 5/8/16 h），**29 nt 原始 + 30 nt（按宿主 K 校正，+2/K）**；读框 = 全基因组参考框图，青绿 = 主翻译框（0），橙红 = +1，紫 = +2；窗口覆盖候选的全部起始。黑线 = PRICE ORF（取最小 P 的 CHX 库）。

### Sup 18　LLC-PK1 与 Vero 的 PEDV P-site 对照
- P`PK1_VERO_ORF_VALIDATION\scripts\corr_figs_editable.py`（A）、`correlation_main.py`（分析：CPM29 = P-site / 同库同读长 PEDV(+) P-site 总数 ×10⁶；**移动块 bootstrap，块长 300 nt，500 次**，给出逐位点与 30-nt bin 的 Spearman 区间）、`nx_local_16h_v11_editable.py`（B：N-x1、N-x2、ORF1b-x13、S-x1 的 Vero/PK1 局部对照；只有计数、份额、比值，**不计算 P 值**）。PK1 底层：P`PK1_psite\scripts\`（由 `riboseq_GIb_mapping_psite` 改写，**offset 在猪宿主 Sscrofa11.1 上独立重算**，Vero 冻结值 29/30 nt = 12 作对照列）。仅含 16 h 两对：PK1 CHX 16 h vs Vero CHX 16 h、PK1 HAR 16 h vs Vero HAR 16 h。

### Sup 19　N-x1 PRM
- P`PRM_Nx1_redraw\scripts\`（`01_calibrate` → `03_digitize` → `04_seq_checks` → `05_redraw`）。**色谱曲线是从厂商图数字化得到的近似，不是原始强度**（脚本 docstring 明说）；C 面板为 N-x1（112 aa）上 PRM signature peptide **LNPLTWVSER（aa 87–96）** 的位置与序列核对（`04_seq_checks.py`）。

---

## 5. 没有（或没有完整）代码的部分

1. **Fig 1A**：示意图（BioRender 类）。
2. **Sup 1**：Western 胶图 + qPCR（`Sup Figure 1.pzfx`）。
3. **Fig 6G**：`fig6.pzfx`（Prism）。
4. **Fig 5 的 GPT 脚本**：`plot_Fig5_stylematched_*.py` 由 GPT 提供，原位置 `<local_downloads>\PEDV_Fig5_plotting_scripts_for_Claude_verification_v2`；本项目已收入 FR`Figure_5\scripts\gpt_plotting_v2\`，并修复过 "`RNA_cov_ORF_body` 缺列时回退第一列" 的 bug（v2）；`Fig5_recolor_Fig4palette` 在其上改配色。
5. **Sup 5**：`plot_canonical_variants.py` 只在数据包里。
6. **Sup 14**：FRASynPlot 外部流程。
7. **Ribo 预处理（cutadapt/bowtie2）与 RNA-seq QC trimming（`analysis_sgRNA_junction/scripts/run_qc_align.sh`）**：只在服务器上，本地没有脚本，只有汇总。
8. **Illustrator 手工步骤**（改标签、删面板、拼版）：不在代码里。

---

## 6. 发现的问题与待核对项（按重要性）

1. **Fig 6A–D 终稿与脚本输出不一致**：C 含 ORF3（代码排除）；B 文字 "max :" 数字缺失；A/C 的卡方 p 值没标。需要以终稿为准重新对账源数据，避免图和源数据对不上。
2. **Fig 3 数值差异与面板来源**：D 的微同源增量与 `Fig3_statistics_methods.md` 略有不同；8 面板拼图找不到单一生成脚本。
3. **命名冲突**：Sup 8/9 的 "N-x4、M-x2、M-x3" 是非经典 sgRNA 衍生 ORF，与 Fig 6 的 +1 框 family "N-x4" 重名但不同物。
4. **Sup 4、9、13 的映射只是"推定"**：图上标签（In vivo / LLC-PK1 等）在 Illustrator 里改过，脚本里没有对应字符串。
5. **两套谱系分型**：Sup 16（S 基因投影，GI-a 23 / GI-b 30）与早期 `genotype_normalized.tsv`（GIa=0、GIb=168）互相不一致，图里只能用前者。
6. **读长集合与 CHX 16 h 样本在不同图里不同**：Fig 4/Sup 10–13/Sup 17 = host_criteria + 813；Fig 6E–F = PRICE/HARTIS 29 nt；Sup 15 的 Ribo 密度用 `Fig6_canonical_background_audit.tsv`（strict_core 29/30 nt，早期 CHX 16 h 为 065）；**方法部分需要逐图写明**。
7. **Sup 14 "Mnvelope" 拼写错误**。
8. **FR 数据包过时**：`final_figure_reproducibility` 只含旧编号的 17 张图（Figure 1–6 + Sup 1–10），没有当前的 Sup 11–19；`00_sample_metadata` 里 `in_current_working_set` 列也过时。

---

## 7. 脚本位置速查

| 主题 | 目录 |
|---|---|
| RNA-seq 比对、junction 提取与分类 | P`GIb_rnaseq_reanalysis\scripts\` |
| Fig 1 MOI 重绘 | P`Fig1_MOI_regrouped\scripts\` |
| Fig 3 | P`GIb_rnaseq_reanalysis\figures\Fig3_work\` |
| Sup 3、6、7（junction 动态） | P`Fig2_Supp_junction_dynamics\scripts\` |
| Sup 4、9（跨宿主） | P`public_data_replication\ESV\{biological_validation_v1,gII_cross_strain_validation_v1}\scripts\` |
| Sup 8（非经典 ORF 类） | P`Noncanonical_reclass\`、P`Noncanonical_ORF_classes\` |
| Ribo 比对 / P-site / 读长集合 | P`riboseq_GIb_mapping_psite\scripts\`、P`module1_rpf_qc\scripts\` |
| Fig 4、Sup 10–12 | P`fig4_v2_host_criteria\`、P`read_length_frame_periodicity\`、P`shared_riboseq_metrics\`、P`supfig12_v2_ABC\`、P`Supplementary_Figure_13\code\` |
| Fig 5、Sup 13 | P`fig5_terminal_ORF\scripts\`、P`fig5_gRNA5prime_figure_package\`、P`Fig5_recolor_Fig4palette\`、P`Fig5_3UTR_translation_audit\` |
| ORF calling（PRICE / Ribo-TISH 等） | P`orf_calling_GIb\`、P`orf_calling_GIb_813\` |
| HARTIS | P`HARTIS_v3\` |
| Fig 6A–D、Sup 15 | K`Result6_genomics_figures\scripts\` |
| Sup 16 | K`Supplementary_Figure_L\scripts\` |
| 序列库与早期比较基因组 | K`N1_conservation_analysis\`、`<PEDV_project>\PEDV_ribo\` |
| Sup 17 | P`Supplementary_Figure_N\PRICE29_locus\code\` |
| Sup 18 | P`PK1_VERO_ORF_VALIDATION\`、P`PK1_psite\` |
| Sup 19 | P`PRM_Nx1_redraw\` |
