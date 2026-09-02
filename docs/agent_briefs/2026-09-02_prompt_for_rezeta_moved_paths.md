# Prompt for the rezeta session: upstream paths moved (2026-09-02)

Hand the block below to the `~/projects/rezeta` session verbatim. It is self-contained.
rezeta vendors nothing from upstream, so this is much lighter than the SFT-WL-B brief:
one live code path (which still works), one broken documentation path, and one warning.

---

上游论文项目 `~/Documents/MyDrafts/STF_lensing/` 在 2026-09-02 做了一次可复现性重组。
你那边引用它的地方很少，而且**没有任何东西损坏**，但有一处文档路径失效、一处值得注意的陷阱。
上游仍然是**只读**的，不要往那边写任何东西。

## 结论先说

**你现在的代码不需要改就能继续跑。** `analysis_emit_vertex_table.py` 读的那个 config 路径
仍然有效，我用你的 `production_cosines()` 逻辑逐字跑过一遍，仍然解析出 40 个分离角
（0.5005′ 到 5000′），与之前一致。

## 权威文件

`~/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/MOVED_PATHS.md`
（旧路径到新路径的完整对照表）。配套：同目录的 `REPRODUCE.md`、`CHANGELOG_REORG_2026-09.md`。

两条结构性变化：`driver_field_emulators/` 已整体并入 `SFT-lensing-paper-analyses/`
（那是个 git 仓库，已推送到 GitHub）；`SFT-lensing-paper-analyses/_archive/` 已退役到废纸篓。

## 需要你改的地方

### 1. 一处失效的文档路径

`docs/history/handoff_prompt_driving_field_channels.md:473`：

```
旧：~/Documents/MyDrafts/STF_lensing/driver_field_emulators/notes/input_ready_b_to_zeta_spec.md
新：~/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses/docs/fk_audit_2026-08/notes/input_ready_b_to_zeta_spec.md
```

同一行附近的 `sections/appendix.tex`（第 471 行）**路径和行号都没变**，正文未被改动。
第 441 到 442 行提到"包括 `driver_field_emulators/` 和 `SFT-lensing-paper-analyses/`"，
现在只剩后者；只读约束不变。

`docs/history/` 下的文件是历史记录，如果你倾向于保持原貌，加一行"2026-09-02 后该路径已变更"
即可，不必逐条改写。

### 2. `SHARED_PROGRAMME.md:18` 与 `findings_S1_sign_chain.md:56`

`SHARED_PROGRAMME.md` 只说"第四个位置是上游论文分析"，仍然成立，可以补一句入口是
`SFT-lensing-paper-analyses/REPRODUCE.md`。

`findings_S1_sign_chain.md:56` 指的是 SFT-WL-B 的 mirror 副本
（`SFT-WL-B/vendor/stf_lensing_mirror/driver_field_emulators/code/callable_fixed/`），
那是他们仓库里的自包含冻结副本，**没有移动，仍然有效**。不要因为上游改名就去改它。

## 一个值得注意的陷阱

`analysis_emit_vertex_table.py` 里的 `CONFIG` 指向

```
.../sftwick_outputs/2PCF/C_corr_op_K_limber_FK/config_L2.yaml
```

这份 config 现在标注为 **SUPERSEDED**：它是六月那次 ell_max=1000、顶点被冻结在角点的运行。
论文实际引用的 FK 来自新的运行目录

```
.../sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/config_L2.yaml
```

**对你当前的用法这没有影响**：你只从 config 里取 40 个分离角的几何，而这个几何在所有运行之间
逐字节相同（按构造如此，我核对过两份 config 解析出的 `positions_grid` 完全一致）。所以继续读旧
的那份是安全的。

但如果你哪天改成读新的那份，注意我在重组时一度让它由 `yaml.safe_dump` 写出，
`positions_grid` 变成了块状序列，你那条内联式正则会**静默匹配到 0 行而不是报错**。
我已经把它改回和其它三份一致的内联写法并推送，现在四份 config 你的正则都能解析出 40 行。
如果你的机器上是旧副本，`git pull` 一下上游包即可。

建议顺手把 `production_cosines()` 加一句断言，例如 `assert len(rows) == 40`，
这样任何格式变化都会立刻暴露而不是悄悄产生空结果。

## 数字上的注意事项

上游同一个量（FK 对 xi_kappa 在 0.5′ 的贡献）在不同文档里有**五个不同的数值**，都不是笔误，
而是四次修正的不同阶段：`+3.086e-5`、`+1.219e-4`、`+4.29e-6`、`+3.84e-6`；论文用的是
`+1.9480e-5`（ell_max=15360，Order-0 的 2.31%）。对照表在

```
SFT-lensing-paper-analyses/docs/fk_audit_2026-08/FK_BASELINE_NUMBERS.md
```

引用 FK 幅度时**必须带上 ell_max**：仅截断本身就差 4.5 倍。上游 13 篇笔记已加了带日期的更正块。
你的 `input_ready_b_to_zeta_spec.md` 引用（上面第 1 条）正是其中之一，它里面的
`+4.29e-6` 是 cut1000 的值。

## 不要做的事

* 不要往 `~/Documents/MyDrafts/STF_lensing/` 写任何东西，整棵树只读。
* 不要 fan out canoes 的 kappa3 HIGH builds（上游因此有过两次 OOM 重启）。
* 不要因为路径变了就去改 SFT-WL-B 的 mirror 副本，那是他们的冻结副本。

## 验证

改完跑一遍你自己的 gate，结果不应有任何变化。如果想确认上游没坏：

```bash
cd ~/Documents/MyDrafts/STF_lensing/SFT-lensing-paper-analyses
python reproduce/check_figures.py     # 17 张图应全部 IDENTICAL
python reproduce/check_numbers.py     # 35 个引用数字对照其产物
```
