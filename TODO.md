# birt 待辦（2026-10-07）

姊妹套件 birtRcpp 的待辦在 `../birtRcpp/TODO.md`。

## 定位

- **主菜：ELGM。** `birt(model, data)` 預設用 RTMB 以 ELGM 配適（`engine = "rtmb"`，`rtmb_control(method = "elgm")`）：
  - 受試者以 AGHQ 積掉，其餘參數用外層的適應性積分；
  - 使用 `dpriors()` 的先驗，並給出邊際概似。
- **同一引擎的最大概似：** `rtmb_control(method = "aghq")`，提供逐人 score，給 `score_test()`、`quantile_score()`、`dif_tree()` 用。
- **JAGS（選用）：** `engine = "jags"`，用於 spike-and-slab 先驗與 ALD 分位數；R2jags 放在 Suggests。
- **核心功能：**
  - `birt()`、`rtirt_syntax()`；
  - 有序題 `ordered =`／`itemtype =`；
  - MNLFA 的 `E()`／`V()`；
  - `hgrm()`（只用 RTMB）；
  - `quantiles()`／`quantile_test()`／`quantile_score()`／`quantile_information()`；
  - `cond_reliability()`；（`rt_information()` 與 `censor =` 屬非核心，見 2b）；
  - `score_test()`、`dif_tree()`；
  - `check_sbc()`、`equations()`、`algorithm()`。
- **不做：** 自寫抽樣器與 ECM（屬於 birtRcpp）。
- **保留：** 和 birtRcpp 共用泛型的轉接層。

## 1. ELGM 主菜要補的

1. [x] **教學改成以 ELGM 為主線**（2026-10-06：教學、範例 ex01–ex10、vignette）。
   - 目前第 1–7 節是 JAGS，用 `fast`／`hopt` 明寫 `engine = "jags"`。
   - 改寫順序：ELGM（後驗、邊際概似、`compare()`）→ 最大概似（score 工具）→ JAGS（ssp、ALD）。
   - 每段輸出都要逐行解釋。
2. [x] vignette 與範例 ex01–ex10 改成以 ELGM 為主（2026-10-06）。
3. [x] ELGM 的速度與穩定性：
   - 預設改成 ELGM 之後，大型模型（多題、雙因素、SHASH）的時間要量測；
   - 外層積分的維度（`hyper = "auto"`）是否需要上限或提示。
   - 2026-10-07：N 1000–2000、20–40 題、81 個外層節點時 ELGM 是 ML 的 1.4–2.0 倍（如 40 題 N = 2000：ML 313 秒、ELGM 488 秒）；`auto` 在 SHASH speed 加 skew 與 scale 調節時選了 4 個方向（625 節點，ELGM 是 ML 的 6 倍）。改成 `auto` 最多 3 個方向並印出訊息：3 個方向的後驗平均數差 ≤ 0.034 SD、SD 差 ≤ 2.1%、logML 差 0.02，外層時間 1/5。`results/birt_groupC/82_timing.csv`、`84_elgm_hyper_dims.txt`。
4. [x] ELGM 版本的 score 工具：目前 `score_test()`／`quantile_score()`／`dif_tree()` 需要最大概似的配適。評估能否在後驗眾數上使用，或維持「先用 aghq 再檢定」的流程，並在文件寫清楚。
   - 2026-10-07：維持「用 aghq 配適再檢定」：眾數上的 score 加總為負的 log 先驗梯度，MAP 版本要置中且型一錯誤率依共變異數估法而變（Debelak et al., 2022）。錯誤訊息與說明文件寫明重配與貝氏做法（有無效果的 ELGM 用 `compare()` 算 BF）；示範：y1 有 DIF 時 BF 2.9e9，y2 無 DIF 時 BF 5.8 支持無 DIF（`83_elgm_bf_dif.txt`）。
5. [x] `hgrm()` 預設跟著 `birt()` 改成 ELGM，但 DIF 流程（`score_test()`）要用 `control = rtmb_control(method = "aghq")`。教學 8.6 要說明這一點。（2026-10-06：8.7 已用 aghq；8.6 加了一句說明。）

## 2. 功能

1. [x] `cond_reliability()`、`rt_information()` 支援有序題與調節載荷（目前會拒絕）。（2026-10-07：在 `at` 的調節變項值計算；與 `ordinal_information()` 一致到 1e-8，與數值 Fisher 資訊一致到 1e-6。）
2. [x] `jags_code()`／`rtmb_code()` 從語法產生程式碼時，接受 `ordered =`／`itemtype =`。（2026-10-07：`rtmb_code()` 也寫出有序題、E() 調節與設限；GRM／GPCM／TPPCM 的腳本重現 −2logL 到 1e-5。）
3. [ ] **延後（非核心，2026-10-05 決定）：調節版 GPCM／TPPCM。**
   - 已經做好的：程式碼、未調節版和 mirt 的比對、暴力積分。
   - 還沒做的：調節版的 SBC、參數回復、TPPCM 步驟鑑別度的調節。

## 2b. 非核心功能（2026-10-07 決定：保留程式碼，不當主打功能）

1. [x] **作答時間的設限概似**（時間上限、把快速作答只當成「很快」）：RTMB 加上 censored likelihood，讓 `rt_information()` 的設計結論可以直接配適驗證。（2026-10-07：`birt(..., censor = list(t1 = c(upper = 60)))`，ELGM 與 ML、常態與 SHASH；JAGS 拒絕。−logL 與暴力積分一致（常態 1e-13）；模擬 N = 1000、18% 設限：設限概似無偏，直接當觀察值的截距偏 −0.066、speed SD 偏 −0.079；`cond_reliability()` 在設限配適上給出限制下的資訊，`rt_information()` 的 `limits` 與設限概似的曲率一致到 1e-5。`results/birt_groupC/81_censor_sim.txt`。）
2. `rt_information()` 的 RT 分位數帶／`limits` 分解也屬非核心（保留）。

## 3. JAGS 側（選用引擎）

1. [ ] 調節版 GRM（ssp）的 SBC 不收斂（`results/birt_refactor/sbc/sbc_summary.md`）：R-hat 中位數 1.2–1.8。
   - 先試 non-centred 改寫（`f = beta*z + sd*u`），讓 beta 和人參數在同一個 glm 區塊更新。
   - JAGS module 的評估（2026-10-05）：成本高，和 birtRcpp 的定位重疊，先不做。
2. [ ] 預設階層載荷先驗（mu ~ N(0, var 10)）的先驗預測會出現 > 25 的載荷：用 prior predictive 重新評估。ELGM 也用同一個先驗，所以這項也關係到主菜。
3. [ ] SBC B 第 71 次複本有一條鏈停住超過 20 分鐘，原因未查。

## 3b. 改寫教學時發現的問題（2026-10-06）

已修：
- `rtmb_code()` 在相關的雙因素模型少了括號，產生的程式無法重現 −2logL（加了回歸測試）。
- `posterior_draws()` 在 `readRDS()` 讀回的配適上回傳空矩陣：恢復 `importFrom(coda)`，讓 coda 隨 birt 載入。

待修：
1. [x] 有序題 ssp 的 `$selection`：`prior_incl` 顯示 .500，但 π 是學出來的（位移約 .40、載荷約 .90），所以 BF10 用的先驗勝算是 1。
   - 2026-10-06：估計表帶上 ssp 組別（`load`、`reg`、`dlt_z`、`alm_z`），`prior_incl` 與 BF10 都用該組 π 的後驗平均數（BF10 原本已用；只有 `prior_incl` 顯示錯）；ex06 印出 .372／.876；測試在 test-hgrm（JAGS）。
2. [x] `algorithm()` 對 ELGM 標「−2 log L」，但那個值含有 −2 log prior；`convergence()` 有寫明，`algorithm()` 沒有。（2026-10-06：標成「-2 log L - 2 log prior」，測試核對數值 = −2·logpost。）
3. [x] 外層質量的警告在 `k_hyper = 3` 時永遠不會出現（edge0 > .49）。
   - 2026-10-06：改成逐方向：最外兩節點的質量 > max(寬 1.5 倍的常態後驗在該規則下的質量, .025)（門檻為啟發式，見 `R/elgm.R`）；常態後驗寬 s 倍時，k = 3、5 在 s > 1.5、k = 9 在 1.76、k = 11 在 2.09 發出警告。測試：合成網格（k = 3–11，1、2、5 維）、N = 200 不警告、N = 25（lsd2 尾巴長，k = 5 時 16%）警告。教學 6.1 的 SHASH 例子：k = 9 時 eps2 3.4% 警告，k = 13 時 1.4% 不警告。
4. [x] 顯示問題（2026-10-06：段名只在適用時加註；`evidence` 用 include.lowest，Inf 為 strong；`plot()` 訊息改成「location shift」；各有測試）：
   - GRM 的摘要印出「Thresholds (partial credit …)」和「Distributions (SHASH …)」；
   - BF10 = Inf 時 `evidence` 欄是空的；
   - 同質變異的 SHASH speed，`plot()` 寫成「normal residuals」。
5. [x] `dif_tree()` 的分割變項是數值時，會窮舉所有切點：2PL 要 3.5 分鐘，N = 500 的聯合模型超過 10 分鐘。考慮預設先分組，或用 `mob_control()` 的選項限制切點。
   - 2026-10-06：新參數 `nbins = 10`：數值變項切成分位數組（有序因子），用 `mob_control(ordinal = "L2")`（maxLM_o；Merkle, Fan & Zeileis, 2014）；`nbins = NULL` 保留窮舉。2PL、N = 600、真切點 x = 6：26 秒，切在 6.01（窮舉：926 秒，切在 6.05；兩次執行時電腦負載不同，窮舉那次較忙）。測試在 test-dif-tree。
6. [x] 教學中 ALD 變異數的式子 σ²(1 − 2p + 2p²)/(p²(1 − p)²)（p = .1 時約 101σ²）要對照 birt 的參數化再確認。（2026-10-06：式子正確，Yu & Zhang, 2005；模擬 2e6 次 101.3 對 101.2。教學改寫：birt 以殘差 SD 為參數、σ 由它推得，所以各 p 的變異數相同；極端分位數不穩的理由改為 p(1 − p)/f(q_p)²（Koenker, 2005）。）

## 4. 文件與整理

1. [x] 教學補上 `quantile_test()`、`quantile_score()`、`rt_information()`、`dif_tree()`、`ordinal_information()` 的小節。（2026-10-06 改寫時已完成：6.4、7.2、6.6、7.3、8.4。）
2. [x] 範例 ex06（hgrm）與 ex10（有序 MNLFA）合併。（2026-10-06：`ex06_ordinal_mnlfa.R`，ELGM → aghq 的 score_test → ELGM → hgrm() → JAGS ssp；範例改為 ex01–ex09；Rscript 跑完。）
3. [x] `tests/testthat/test-review.R` 的回歸測試併入各主題測試檔。（2026-10-06：分到 test-syntax、test-rtmb、test-hgrm、test-fit。）
4. [x] NEWS 改寫成首次發布的功能說明。（2026-10-06）

## 5. 研究線（以 birt 的工具為主）

### 5a. 分位數（2026-10-06 盤點；文獻與驗證：`../notes/qr_decision_indices_literature.md`）
| # | 題目 | 功能 | 現況 |
|---|---|---|---|
| 1 | speed 的潛在分位數迴歸：ALD 對 SHASH（論文第 4 章） | `ald()`（JAGS）、`shash()` + `quantiles()`（ELGM／ML）；對照 birtRcpp `rtirt_latent(quantile =)` | ALD 在 p = .25 的偏誤來自潛在 ALD，待解 |
| 2 | 分位數 cross 關係（論文第 5 章） | `quantile_test()`、sup-t 同時信賴帶；對照 birtRcpp `rtirt_cross(quantile =)` | 可以跑，已驗證 |
| 3 | 測量誤差下的分位數迴歸 | `quantiles()` 對照 `scores()` 加 quantreg | 需要寫模擬 |
| 4 | 分位數效果的不變性 | `score_test(quantile_score(fit, p), group)` | 已驗證 |
| 5 | RT 的訊息量落在哪裡、時間上限怎麼設 | `rt_information()`、`cond_reliability()` | 已驗證；需要實際資料 |
| 6 | 快速猜題與快尾 | `rt_information(limits = list(fast =))` | 設限概似已有（2b.1，`censor =`） |
| 7 | SHASH 設錯時的穩健推論 | `quantile_score()$se.robust` | 需要寫模擬 |
| 9 | 哪一條分位數迴歸對 ability 的訊息最大（cross 關係 `t ~ speed + ability`） | `quantile_information()`：I_τ = f(q_τ)² β_a(τ)² / (τ(1 − τ))（Godambe 資訊），去掉 speed，對整個 RT 分配的效率 | 2026-10-07 完成；常態閉式解、SHASH 分位數函數（1e-5）、二分 RT 的 Monte Carlo（3% 內）；設計 `../notes/quantile_information_design.md` |
| 8 | 分位數效果的貝氏推論：ELGM、JAGS ALD、birtRcpp Gibbs ALD | `quantiles()` 在 ELGM 配適上 | 可以跑 |

- 建議：第 2 與第 5 項合寫。
- 不做：「分位數迴歸係數的信度」（2026-10-06 決定）。

### 5b. 其他
- RT-IRT 信度論文：等使用者提供資料（`../PLAN.md`）。
- 分數檢定的保守時間順序版本：未定。
