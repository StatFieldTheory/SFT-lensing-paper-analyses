# Prompt for the SFT-WL-B session: upstream paths moved (2026-09-02)

Hand the block below to the `~/projects/SFT-WL-B` session verbatim. It is self-contained.

---

上游论文项目 `~/Documents/MyDrafts/STF_lensing/` 在 2026-09-02 做了一次可复现性重组，你引用和 vendor 的路径大部分变了。请据此更新本仓库的引用。上游仍然是**只读**的，不要往那边写任何东西。

## 权威文件

`~/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/MOVED_PATHS.md`

旧路径到新路径的完整对照表，并逐一标注了你 vendor 的每个文件字节是否变化。先读它，再动手。配套还有同目录的 `REPRODUCE.md`（怎么复现每张图和每个数字）和 `CHANGELOG_REORG_2026-09.md`（做了什么、怎么回滚）。

## 两条结构性变化

1. **`driver_field_emulators/` 不存在了。** 它是 2026-08 的工作树，现已整体并入
   `SFT-lensing-paper-analyses/`（那是个 git 仓库，已推送到
   `github.com:StatFieldTheory/SFT-lensing-paper-analyses`）。生产 FK 路径现在就在复现包里。
2. **`SFT-lensing-paper-analyses/_archive/` 也不存在了。** 2026-05 的 canoes pipeline 与
   2026-06 的清理归档已退役到废纸篓（附 md5 清单）。其中仍被使用的三样东西已 vendor 进活动树。

## 需要你改的地方

### A. `vendor/*/PROVENANCE.md` 的 upstream 路径列

你 vendor 的 14 个上游文件里，**11 个字节未变**，钉住的 SHA-256 仍然成立，只需要改路径列：

| 你的 vendor 名 | 新的上游路径 |
|---|---|
| `perm_aware_kappa3_callable.py` | `SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/equal_time_limber_cut15360_permaware/perm_aware_kappa3_callable.py` |
| `callable_fixed/README.md` | 同上目录 `README.md` |
| `callable_fixed/table_permclosed.npz` | 同上目录 `table_permclosed.npz`（原 `products/table_permclosed_cut15360.npz`） |
| `equal_time_limber_kappa3_callable.py` | 路径未变 |
| `equal_time_limber_kappa3_z5covgrid16_...npz` | 路径未变 |
| `xi_C_corr_op_O0.npz`、`xi_C_corr_op_K_limber_FF.npz`、`xi_C_corr_op_K_limber_FK.npz` | 路径均未变 |
| `plot_cl_EB_polarization.py`、`_plot_style.py` | 路径均未变 |
| `cutoff_study/cutoff_vertex.npz` | `SFT-lensing-paper-analyses/sachs_sft/analyses/fk_audit_2026-08/products/cutoff_vertex.npz` |

**3 个字节变了**，如果你重新 vendor，SHA 会对不上：

| 文件 | 旧 SHA-256(16) | 新 SHA-256(16) | 变了什么 |
|---|---|---|---|
| `compare_callables.py` | `73ac4f587532924d` | `b9f7989e0c0b7754` | 硬编码的 Order-0 路径改成相对路径 |
| `tests/test_perm_aware.py` | `1eaf330b5712081b` | `d0cbaaa51ee05a89` | 定位表文件的 `parents[...]` 层数 |
| `plot_analysis3_cl_decomposition.py` | `455dd479c0278147` | `f28f9f5d21424bc9` | **语义变化，见 B** |

顺带一提：`table_permclosed.npz` 现在在上游就叫这个名字（即 callable 自己解析的名字），所以你
mirror 里那步重命名不再需要了。

`vendor/consumer_pk/PCAMBz0.txt` 来自 canoes 而非本项目，不受影响。
`vendor/bihalofit_ref/` 来自 GitHub，不受影响。

### B. 一处语义变化，别只当成改路径

`plot_analysis3_cl_decomposition.py` 的 `FK_NPZ` 默认值变了：

```
旧：sftwick_outputs/2PCF/C_corr_op_K_limber_FK/xi_C_corr_op_K_limber_FK.npz
新：sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz
```

旧的是六月那份 ell_max=1000、顶点被冻结在角点的 sweep；新的是论文实际引用的收敛结果，
0.5′ 处大 4.5 倍。**你 vendor 的数据文件本身没变**，所以现有 gate 的结果不受影响；但如果你
重新 vendor 这个脚本并依赖它的模块级常量，读到的 FK 就换了一份。旧 sweep 仍在原地保留（其
README 加了 SUPERSEDED 标注），显式指名它的代码继续可用。

用新 sweep 重算出来的 B/E 就是论文引用的那组：ell=60 处 0.95，降到 ell=1500 处 0.42，
50≤ell≤1500 区间中位数 0.629。

### C. `docs/decision_log_bmode.md`

* **第 15 行**：`[R]` 这个引用键定义为 `~/Documents/MyDrafts/STF_lensing/driver_field_emulators/`，
  该目录已不存在。请重新定义或删除这个键，并检查所有用到 `[R]` 的条目。
* **第 206 行**：写着生产 callable 的 `np.sort(...)` 加 `np.round(..., 8)`
  "Still unfixed in production"。**这句现在不成立了**：那个 callable 已不是生产路径。
  生产顶点是 `equal_time_limber_cut15360_permaware/`（permutation-aware，ell_max=15360），
  旧文件旁边放了 `DEPRECATED.md` 说明。文件本身和行号都没变，只是"生产"的含义变了。
* **第 2367 行**：`driver_field_emulators/notes/finding_grid_artifact_measured.md`
  → `SFT-lensing-paper-analyses/docs/fk_audit_2026-08/notes/finding_grid_artifact_measured.md`
  （`notes/` 下所有文件都按原名迁到这个目录）。
* **第 2376 行**指的那份 FK sweep 仍在原路径，但已标注 SUPERSEDED。
* 所有 `[P] sections/...` 与 `[P] .../mathematica/...` 的引用**行号和路径都没变**，正文未被改动。

### D. 其它引用了旧路径的文件

`analysis/wpb1/b11_measure_exponent.py`、`analysis/wpb1/FINDINGS_B11.md`、
`analysis/wp0/extract_collapsed.py`、`docs/handoff_canoes_bmod_phase.md`、
`docs/handoff_rezeta_full_review.md`、`docs/plan_bmode_programme.md` 里都提到了
`driver_field_emulators/...`。常见的两个对应关系：

* `driver_field_emulators/code/rebuild/assemble.py`
  → `SFT-lensing-paper-analyses/sachs_sft/callables/kappa3_vertex/rebuild/assemble.py`
* `driver_field_emulators/notes/input_ready_b_to_zeta_spec.md`
  → `SFT-lensing-paper-analyses/docs/fk_audit_2026-08/notes/input_ready_b_to_zeta_spec.md`

`docs/handoff_stf_lensing_cleanup.md` 是你当初写给上游的那份委托书，它描述的是**重组前**的状态。
建议保留原文并在顶部加一行说明"已于 2026-09-02 执行完毕，现状见上游 CHANGELOG_REORG_2026-09.md"，
不要逐条改写它。

### E. 一个数字上的注意事项

上游同一个量（FK 对 xi_kappa 在 0.5′ 的贡献）在不同文档里有**五个不同的数值**，都不是笔误，
而是四次修正的不同阶段：`+3.086e-5`、`+1.219e-4`、`+4.29e-6`、`+3.84e-6`，论文用的是
`+1.9480e-5`（ell_max=15360，Order-0 的 2.31%）。对照表在
`SFT-lensing-paper-analyses/docs/fk_audit_2026-08/FK_BASELINE_NUMBERS.md`。
引用 FK 幅度时**必须带上 ell_max**：光是截断本身就差 4.5 倍。上游有 13 篇笔记已加了带日期的更正块。

## 不要做的事

* 不要因为路径变了就重新 vendor 全部文件。11 个字节未变，钉住的哈希仍然成立；重新 vendor 只会
  让"逐字节一致"这个已验证的声明失去意义。真要重新 vendor，就重新验证并更新哈希。
* 不要动 `vendor/stf_lensing_mirror/`（`chmod 555`、带封闭完整性声明）的内容。它是自包含副本，
  仍然照常工作；只有 `PROVENANCE.md` 里的 upstream 路径列过时了。
* 不要试图为 `cutoff_vertex.npz` 找生产脚本。我在整个复现包和已退役的归档里搜过它的数组名
  （`tree_curves`、`tree_edges`、`tree_chi_mpc` 等）与文件名，**确认不存在**。你
  `vendor/cutoff_study/PROVENANCE.md` 里那条警告依然成立：它是个 stamp，只能当形状交叉检验用，
  不能单独承载数值结论。

## 验证

改完后跑一遍你自己的 gate，确认数值没变（本次重组不应改变你的任何结果）。如果想确认上游没坏，
在上游目录下跑：

```bash
cd ~/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses
python reproduce/check_figures.py     # 17 张图应全部 IDENTICAL
python reproduce/check_numbers.py     # 35 个引用数字对照其产物
```

上游三个仓库都已推送并与远端同步，回滚点是各仓库的 `pre-reorg-2026-09` 标签
（论文仓库在 Overleaf 上不支持标签，回滚点是 commit `1f0ec58`）。
