# birt

**BIvariate Response Time Modeling with 'lavaan'-Style Syntax**

用 lavaan 風格的語法寫 IRT／作答時間模型，最多容納兩個潛在變項（例如能力 `ability` 與速度 `speed`）。

**主菜是 ELGM**：`birt(model, data)` 預設用 RTMB 把受試者的潛在變項 θ 以適應性 Gauss–Hermite 積分（AGHQ）積掉，其餘參數再用外層的適應性積分得到近似後驗（extended latent Gaussian model；Stringer, Brown & Stafford, 2023），使用 `dpriors()` 的先驗，並給出比較模型用的邊際概似。

| | 條件概似 p(y \| θ, ψ)：θ 和參數一起抽樣 | 邊際概似 ∫ p(y \| θ, ψ) p(θ \| ψ) dθ：θ 積掉 |
|---|---|---|
| 貝氏 | `engine = "jags"`（選用；spike-and-slab、ALD） | **預設**：`engine = "rtmb"`，`rtmb_control(method = "elgm")` |
| 最大概似 | 不提供（聯合最大概似不一致） | `rtmb_control(method = "aghq")`（或 `"laplace"`） |

- **ELGM（預設）：** 近似貝氏；SHASH 殘差、尺度調節（`V(y) ~ z`）與偏態調節（`shash(skew = )`）建立分配模型，一次配適就由 `quantiles(fit, p)` 得到每個分位數的效果；條件信度 `cond_reliability()`、作答時間的訊息 `rt_information()`。
- **最大概似（`rtmb_control(method = "aghq")`）：** 同一個邊際概似，同 GLMMadaptive、mirt；逐人的 score 給 `score_test()`、`quantile_score()`、`dif_tree()` 使用。
- **JAGS（`engine = "jags"`，選用，需要安裝 JAGS）：** 模型語法翻成一個乾淨的 JAGS 檔（`jags_code(fit)`），提供 spike-and-slab 先驗 `prior("ssp")`、分位數分配 `ald(p)`、邊際概似的 DIC／WAIC／PSIS-LOO 與後驗預測檢查。

三者共用同一份模型語法和 `equations()`（最後一段寫明是哪一種概似）；`algorithm()` 寫出各自的計算。人分數：RTMB 是給定估計值（或後驗）的 EAP，JAGS 是抽樣到的 θ 的後驗平均與 SD。（這裡的「條件」指給定 θ；`cond_reliability()` 的「條件信度」是另一個意思：給定 θ 的信度。）

## 安裝

```r
install.packages("birt")                       # CRAN 上架之後
# install.packages("remotes")
remotes::install_github("jiewenTsai/birt")     # GitHub 上的版本
```

RTMB、lavaan 等必要套件會一起安裝。

JAGS 引擎是選用的：要用 `engine = "jags"` 時，先安裝 [JAGS](https://mcmc-jags.sourceforge.io/)（≥ 4.3）與 `install.packages("R2jags")`。

教學與範例：[birtExamples](https://github.com/jiewenTsai/birtExamples)。套件裡也有教學原始檔 `inst/tutorial/birt_tutorial.qmd`，以及範例腳本 `inst/examples/ex01`–`ex09`（`system.file("examples", package = "birt")`）。

與 birtRcpp 共用存取函數的名稱（`estimates`、`scores`、`reliability`、`fit_indices`、`convergence`）。兩個套件同時載入時，不論誰後載入都能正確派發。

## 快速開始

```r
library(birt)
d <- sim_rtirt(N = 500, K = 10)

model <- rtirt_syntax(paste0("y", 1:10), paste0("t", 1:10))   # 產生語法，也可以手寫
model
#> ability =~ y1 + y2 + y3 + y4 + y5 + y6 + ...
#> speed =~ -1*t1 + -1*t2 + ...
#> ability ~~ speed
fit <- birt(model, d)                       # 預設：RTMB + ELGM（近似貝氏，邊際概似）
summary(fit)                                # 後驗平均、後驗 SD、95% 可信區間、prior 欄
ml <- birt(model, d, control = rtmb_control(method = "aghq"))   # 同一個邊際概似的最大概似

# JAGS（選用）：spike-and-slab 交叉負荷
fj <- birt(rtirt_syntax(paste0("y", 1:10), paste0("t", 1:10), cross = "ssp"), d, engine = "jags")  # 4 條鏈平行、10000 次、燒 5000
```

真實資料常見的寫法（原始秒數、受試者 ID）：

```r
items <- grep("^ME[0-9]+$", names(raw), value = TRUE); times <- paste0(items, "_S")
model <- rtirt_syntax(items, times, cross = "free")
fit <- birt(model, raw, log_rt = TRUE, id = "IDSTUD")          # 先取 log；scores() 帶 ID
```

**`rtirt_syntax()` 的選項：**
- `cross`：`"none"` / `"free"` / `"ssp"`
- `structure`：`"cor"` / `"path"` / `"none"`。有交叉負荷時只能 `"none"`，因為共同的交叉負荷等於相關或路徑。
- `cov` 與 `cov_on`：潛在回歸的共變數，以及要回歸到哪個潛在變項
- `speed_var = "fixed"`：產生 `speed ~~ 1*speed`

**資料防呆：**
- 二元題若是 1/2 這類兩值編碼，會報錯並提示怎麼重新編碼；
- 連續指標看起來是原始秒數（全為正、右偏、中位數 > 10）時，會警告並建議用 `log_rt = TRUE`。

## RTMB 引擎：分配模型與所有分位數

```r
m5 <- paste(rtirt_syntax(items, times, cross = "free"),                # 第 5 章：θ 與各題 log RT 的交叉關係
            sprintf("V(%s) ~ ability", paste(times, collapse = " + ")), sep = "\n")   # log 尺度隨 θ 改變
f5 <- birt(m5, raw, log_rt = TRUE, engine = "rtmb", family = list(.continuous = shash(skew = "ability")), control = rtmb_control(method = "aghq"))
summary(f5)                       # est、SE、z、p、95% CI；殘差分配表（sd、skew、skew~ability、tail）與 V(t) ~ ability
quantiles(f5, p = c(.1, .5, .9))  # θ 對各題 log RT 條件分位數的效果（delta-method SE）
plot(f5, ald = fit_ald)           # 分位數曲線 + 95% 帶；可疊上 JAGS 的 ald(p) 估計

m4 <- paste(rtirt_syntax(items, times, structure = "path", cov = X),   # 第 4 章：speed ~ θ + X
            sprintf("V(speed) ~ %s", paste(c("ability", X), collapse = " + ")), sep = "\n")
f4 <- birt(m4, raw, log_rt = TRUE, engine = "rtmb", family = list(speed = shash()), control = rtmb_control(method = "aghq"))
compare(normal = f4n, shash = f4)                                      # -2logL、AIC、BIC
```

**模型：** 每個潛在變項都是獨立標準常態的轉換。
- 殘差可以是常態或 SHASH（sinh-arcsinh；Jones & Pewsey, 2009），可用在第二個潛在變項（例如 speed）與連續指標上。
- SHASH 經過標準化，使 T(0) = 0、T′(0) = 1，因此截距是條件中位數，`f2 ~ x` 的係數是對中位數的效果。
- 語法 `V(y) ~ z` 讓 log 尺度、`shash(skew = )` 讓偏態成為潛在變項或共變數的線性函數，分位數效果因此可以隨 p 改變。

**估計：**
- 先用 Laplace 近似：從常態模型開始，再估計完整模型。
- 接著做逐人的 AGHQ：二維 k × k 乘積規則，預設 k = 9，放在條件眾數並以 Hessian 縮放。
- 每一輪把節點重新放到新估計值的眾數，再重新評估目標函數；若變差就把步長減半，直到 −2logL 不再改變。
- `rtmb_control(method = "laplace")` 只用 Laplace。

**近似貝氏（ELGM，預設）：** `rtmb_control(method = "elgm")` 把模型當作 extended latent Gaussian model（Stringer, Brown & Stafford, 2023）。

```r
fe <- birt(m, d, engine = "rtmb", dp = dpriors(cross = "dnorm(0, 100)"),
           control = rtmb_control(method = "elgm"))      # 預設：inner = "aghq"，hyper = "auto"
summary(fe)                    # 後驗平均、後驗 SD、95% 可信區間、prior 欄
compare(m1 = fe, m0 = fe0)     # log 邊際概似與 Bayes factor
```

- 高斯場 W 只放人的潛在變項。內層逐人積分：`inner = "aghq"`（預設），或 `"laplace"`（原始 ELGM）。
- 外層對其餘參數的後驗做積分，用混合網格：
  - `hyper` 的每一維放 `k_hyper` 個節點，並排在 Cholesky 分解最前面。預設 `"auto"`：SHASH 潛在變項取形狀參數，否則取潛在相關、路徑與 SD。每一維依 ±2 SD 處的後驗下降量，左右分別調整寬度；
  - 其他方向（題目參數）以 Gaussian 處理，在每個外層節點移到條件後驗眾數（`cond_modes = TRUE`）。
- 先驗：`dpriors()` 與 `prior("...")`，含截斷 `T(,)`；不支援 `prior("ssp")`。
- 報表：後驗平均、SD、95% 可信區間（由積分動差算出；SD 在 log 尺度、相關在 Fisher-z 尺度）。`quantiles()` 與 EAP 分數也都是外層節點上的混合。
- **和 NUTS 的比較**（`scripts/38_elgm_check.R`）：NUTS 用 tmbstan 抽同一個聯合後驗，人與參數一起抽；N = 300–400。
  - van der Linden 模型：ρ 0.431（SD 0.075），NUTS 0.428（0.076）；所有參數的後驗平均都在 NUTS 的 0.25 個 SD 以內。JAGS 引擎（相同 `dpriors()`）與 NUTS 的差距在 0.06 個 SD 以內。
  - 交叉負荷量加 SD 0.1 的先驗：最大差距 0.27 個 SD。
  - SHASH speed：偏態 0.530（0.235），NUTS 0.528（0.247）；尾部 0.982（0.164），NUTS 0.974（0.167）。
  - 限制：外層以 Gaussian 處理的方向，準確度與 Laplace 相同。識別很弱、後驗右尾很長的鑑別度（短測驗）會低估後驗 SD（例如 0.60，NUTS 為 1.02），這時請用 JAGS 引擎。
  - 改用 `inner = "laplace"` 時，最大差距由 0.25 增為 0.45 個 SD。若不做條件眾數、只把 Gaussian 方向放在條件平均線上（pilot 的做法），SHASH 偏態的後驗 SD 只有 0.14（NUTS 為 0.23）。

**分位數效果：** 是條件分位數對各預測變項的導數，在模型隱含的調節變項分配上取平均（AME），標準誤用 delta method。常態且沒有調節時，各個 p 的效果相同，就是均值模型的係數。

**與 JAGS 引擎的差別：**
- 最大概似沒有先驗：`prior()` 修飾會被忽略並提示；`method = "elgm"` 會使用這些先驗（`prior("ssp")` 除外）。
- 不接受 `ald()`：ALD 是單一分位數的 working likelihood，RTMB 引擎改用 `shash()`。
- 小樣本或題數很少時，最大概似可能發散（例如鑑別度 → ∞）。套件會警告，此時請改用 JAGS 引擎。

**TIMSS 驗證**（`scripts/37_birt_rtmb_check.R`）：
- 常態與同質 SHASH 模型的 −2logL，和先前獨立撰寫的 RTMB 腳本完全一致：
  - 第 5 章 normal homo 23470.91、normal het 22932.02；
  - 第 4 章 normal homo 23043.78、SHASH homo 22966.50、SHASH het 22932.47（腳本為 22932.39）。
- 第 5 章 SHASH + θ 偏態：
  - ρ_i(p) 與腳本的相關為 .996，與 birt ALD 的相關為 .950；
  - 速度變異數改為自由時，−2logL 由 22480 降到 21613。
- 第 4 章 θ 的分位數效果（p = .05 至 .95）：.019 至 −.212，腳本為 .020 至 −.214。

## 語法

| 寫法 | 意思 |
|---|---|
| `f =~ y1 + y2` | 負荷量。每個因素的第一行 `=~` 是測量負荷量（正），之後各行是交叉負荷 |
| `f =~ -1*t1` | 固定負荷量 |
| `f =~ r*t1 + r*t2` | 以標籤設定相等 |
| `f =~ prior("ssp")*t1` | spike-and-slab |
| `f =~ prior("dnorm(0, 100)")*t1` | 任意 JAGS 先驗（精確度參數化） |
| `f1 ~~ f2` | 因素相關（預設為 0） |
| `f ~~ 1*f` | 固定（殘差）變異數 |
| `f ~ x1 + prior("ssp")*x2` | 潛在回歸 |
| `f2 ~ f1` | 兩個潛在變項之間的路徑 |

**先驗：**
- `summary()` 在每個係數旁列出它的先驗（`prior` 欄），格式同 blavaan。
- 單一參數用 `prior("...")` 自訂，例如 `prior("dnorm(0, 100)")`、`prior("ssp")`。
- 整類參數的預設用 `birt(..., dp = dpriors(cross = "dnorm(0, 100)"))` 修改，對應 blavaan 的 `dpriors()`。
- 寫法都是 JAGS 的分配，常態為精確度參數化。

**spike-and-slab 的讀法：**
- `p_incl` 是參數不為 0 的後驗機率；`selected` 依中位數機率模型判定（p_incl > .5）。
- `BF10` = 後驗勝算 ÷ 事前勝算。
- 預設的共用納入機率是從資料學的（Beta(1, 1)）。多數效果明顯時它會偏高，邊緣參數也會被選進來；要做檢定請固定：`ssp_p = 0.5`。

**識別規則：**
- 有自由負荷量的因素，變異數固定為 1；負荷量全部固定的因素（例如 speed），變異數自由估計。
- 若 `ability` 對 RT 有自由的交叉負荷，則 `speed ~ ability` 不可識別，套件會拒絕（只對部分 RT 有交叉負荷時仍可識別）。

## 分配（`family`）

```r
family = list(speed = ald(0.25))          # 結構層：speed 的 0.25 分位數回歸
family = list(speed = ald(c(.1, .5, .9))) # 多個分位數，一次跑完，回傳 birt_quantile
family = list(.continuous = ald(0.5))      # 測量層：所有連續指標用中位數模型（穩健）
family = list(lrt_ME62095 = ald(0.5))     # 只改單一指標
```

- **名稱依語法中的角色解析：** 潛在變項名稱對應結構層，指標名稱對應測量層，`.continuous` 代表全部連續指標。二元題固定用 logit。
- **變異數固定的潛在變項**（例如 ability）：尺度由「隱含變異數 = 固定值」決定，所以不同 p 之間 ability 的尺度可以比較。
- **限制：** ALD 不能和 `f1 ~~ f2` 並用，請改寫成 `f2 ~ f1`。

## 有序題與調節敘述（MNLFA）：`ordered =`、`E()`、`V()`

有序題（Likert）用 `ordered = TRUE`（所有整數值的指標）或列出題名。有序題可以和二元題、連續指標（例如 log 作答時間）放在同一個模型，最多兩個潛在變項。調節敘述對所有指標型態都能用。

有序題的模型用 `itemtype =` 選：`"grm"`（graded response，預設）、`"gpcm"`（generalized partial credit）、`"tppcm"`（two-parameter partial credit，每一步有自己的鑑別度）；可以整份一種，或用具名向量逐題指定。三種都能加調節敘述。GPCM 和 mirt 的結果一致；GPCM／TPPCM 目前是實驗性功能（調節版還沒做模擬校準）。飽和的 TPPCM 是 nominal 模型的重新參數化，和 mirt 的 nominal 模型 −2logL 相同。

例：每個人只有一個（整份問卷的）作答時間時，讓時間調節 GRM 的鑑別度與閾值（類似 DIF），以及特質的平均數與變異數：

```r
d$logT <- log(d$Time / 196)                       # log 每題秒數
m <- '
  SA =~ SelAwa1 + SelAwa2 + SelAwa3 + SelAwa4 + SelAwa5 + SelAwa6
  E(SelAwa1 + SelAwa2 + SelAwa3 + SelAwa4 + SelAwa5 + SelAwa6) ~ prior("ssp")*logT      # 閾值（uniform DIF）
  E(SelAwa1 + SelAwa2 + SelAwa3 + SelAwa4 + SelAwa5 + SelAwa6) ~ prior("ssp")*logT:SA   # 鑑別度（nonuniform DIF）
  E(SA) ~ logT                                                                           # 平均數
  V(SA) ~ logT                                                                           # 變異數
'
fit <- birt(m, d, ordered = TRUE, engine = "jags")
summary(fit)          # 題目參數、α / β（p_incl、BF10）、impact、依時間分組的精確度
plot(fit, type = "dif"); plot(fit, type = "information")
```

調節敘述沿用 MNLFA 的寫法：把條件動差對調節變項做迴歸。

- 題目：E(y_i | SA, z) = ν_i(z) + λ_i(z)·SA。類別題（二元、有序）的位移在潛在尺度上：a(SA − b − δ)，δ = βz；連續題是截距位移。
- 潛在變項：E(SA | z) 與 V(SA | z)。
- 連續指標：`V(t1) ~ z` 讓殘差的 log SD 隨 z 改變（z 也可以是潛在變項，限 RTMB）。


| 語法 | 意思 |
|---|---|
| `E(item) ~ z` | 該題閾值隨 z 平行平移（uniform DIF）；連續題是截距 |
| `E(item) ~ z:SA` | 該題鑑別度（載荷）隨 z 改變（log 尺度；nonuniform DIF） |
| `V(t) ~ z` | 連續指標的殘差 log SD |
| `E(item1 + item2) ~ pa*z:SA` | 同一標籤代表效果相等；全部題目共用就是 Ferrando 的受試者精確度模型 |
| `E(SA) ~ z`（或 `SA ~ z`） | 特質平均數（impact） |
| `V(SA) ~ z` | 特質變異數（log SD 尺度） |
| 不寫，或 `0*z` | 定錨題 |

- 修飾詞與其他模型相同：`prior("ssp")`、`prior("...")`、標籤。
- 調節變項照原值使用（要效果是每 1 SD，請先標準化；`hgrm()` 預設會幫你標準化）。
- **識別：** `E(SA) ~ z` 搭配「每一題都有 `E(item) ~ z`」時不可識別（共同平移等於平均數平移），套件會拒絕。至少要留一題定錨題，或用 `prior("ssp")`。

`engine = "rtmb"` 以 AGHQ 做最大概似，可用自由效果和標籤（不支援 `prior("ssp")`），並提供 SE、AIC、BIC；也可用 `rtmb_control(method = "elgm")` 做近似貝氏。對數概似和直接數值積分的差距 < 10⁻⁶。

`hgrm(model, data, moderators, a, b, impact, anchor)` 是捷徑：它把 `a`、`b`（`"free"`、`"common"`、`"none"`）與 `impact` 翻成上面的語法（存在 `fit$model`），再用 RTMB 引擎配適（最大概似，或 `control = rtmb_control(method = "elgm")`）。預設 `a = "common"`、`b = "none"`、`impact = "mean"`。DIF：先用 `score_test()` 逐題檢定，再用 `b = "free"` 加 `anchor` 放寬。spike-and-slab 篩選請直接寫 `birt()` 語法（JAGS）。完整說明見教學文件第 8 節。

## 檢查 JAGS 的翻譯：`check_sbc()`

`check_sbc(model, data)` 做模擬校準（SBC；Talts et al., 2018）：每一次從先驗抽參數、照模型模擬資料（人數、共變數、調節變項、遺漏位置都和 `data` 相同），再用 `birt()` 配適，記下真值在後驗抽樣中的名次。模型、先驗和抽樣器一致時，名次是均勻分配。模擬程式是依模型的定義在 R 裡另外寫的，不是從產生的 JAGS 程式碼來的，所以語法翻成 JAGS 時的錯誤也會被抓到。每一次都是一個 JAGS 配適，正式檢查要 100 次以上。

## Score-based 的參數不變性檢定：`score_test()`

RTMB 引擎的最大概似配適（含 `ordered = TRUE`）可以用 `score_test()` 檢定題目參數是否隨調節變項改變（DIF、題目位置效應），只需要一次配適：

```r
f <- birt(rtirt_syntax(items, times), d, engine = "rtmb", control = rtmb_control(method = "aghq"))
score_test(f, d$gender, by_item = TRUE)          # 性別 DIF，逐題並做 Holm 校正
score_test(f, rowSums(d[times]))                 # 沿總作答時間（聯合模型中對作答參數有效）
```

- 每人的 score（`estfun()`）由 RTMB 對每人的 AGHQ 對數概似取 Jacobian，和逐人數值微分一致到 10⁻¹⁰。
- 調節變項和因素分數有關時（impact），請把它放進潛在迴歸（`f ~ z`，順序題用 `E(f) ~ z`）再檢定；函數會提醒。
- birtRcpp 也有同名函數；兩者都載入時，用 `birt::score_test()` 可以同時處理兩個套件的配適。

## 物件與存取函數

**摘要表物件：** `summary(fit)` 一次算好所有表格，回傳 S3 物件（`birt_summary`）。每張表都是 data.frame，印出時標上名字，用 `s$items` 等取出。版面參考 tppcm 的 `irt_pars()`，參數表用 blavaan 的分段與 Prior 欄。

| 元素 | 內容 |
|---|---|
| `$fit` | 配適指標（DIC／WAIC／LOOIC 與最大 R̂；或 logLik／AIC／BIC） |
| `$items`、`$se` | 每題一列：負荷量、截距、殘差變異數（二元題另列 IRT a、b），以及同版面的後驗 SD 或標準誤 |
| `$parameters` | 每個參數一列，依 lavaan 分段，含 Prior 欄（貝氏）或 z、p（ML）；`as.data.frame(s)` 也會得到它 |
| `$selection`、`$moderation`、`$precision`、`$reliability` | spike-and-slab 選擇、順序題的調節效果、依調節變項的精確度、EAP 信度 |

`print(s, tables = c("fit", "parameters"))` 只印部分表格。

**模型式子：** `equations(fit)` 把模型寫成式子，風格接近 SAS NLMIXED，包含參數符號、適用的題目集合，貝氏配適時另附先驗。預設輸出 Unicode；`equations(fit, "latex")` 輸出 `aligned` 區塊，可以放進論文附錄。

**計算步驟：** `algorithm(fit)` 一步一步寫出配適怎麼算：JAGS 列出完整條件分配與 JAGS 指派的抽樣器（`rjags::list.samplers()`）和 MCMC 設定；RTMB 列出 Laplace／AGHQ 的積分中心、節點、最佳化輪數與 SE；ELGM 列出內外層積分、節點權重、log 邊際概似與後驗摘要的算法。同樣有 Unicode 與 `"latex"` 兩種輸出。

JAGS 的配適結果都屬於 `birt_fit` 類別（`c("birt", "birt_fit")`），共用下表的存取函數。RTMB 的結果（包括 `hgrm()`）是 `birt_rtmb`，支援：

- `summary`、`coef`、`estimates`（est、se、z、pvalue、CI）、`scores`、`reliability`
- `fit_indices`（logLik、Deviance、AIC、BIC）、`convergence`
- `logLik`、`AIC`、`BIC`、`vcov`
- `quantiles`、`plot`（分位數曲線）、`compare`
- 分位數效果的決策：`quantiles()` 的同時信賴帶（`band.lower`、`band.upper`）、`quantile_test()`（效果是否隨分位數改變）、`quantile_score()`（每人的影響函數、穩健 SE、`score_test(quantile_score(fit, p), group)`）
- 信度與訊息：`cond_reliability()`（條件信度）、`rt_information()`（作答時間的訊息落在 RT 分配的哪一段、時間上限保留多少）、`ordinal_information()`（有序題的測驗訊息）
- 不變性：`score_test()`（逐題 DIF，一次配適）、`dif_tree()`（DIF 樹，需要 partykit）

| 函數 | 內容 |
|---|---|
| `summary(fit)` | 分段參數表：負荷量、回歸、變異數、截距、ssp、信度、適配 |
| `estimates(fit)` | lavaan 式參數表（est, sd, q025, q975, rhat, ess, scale, p_incl, BF10） |
| `coef(fit)` | 自由參數的後驗平均數（具名向量） |
| `scores(fit)` | 每個人各潛在變項的 EAP 與後驗 SD（`id =` 給的 ID） |
| `reliability(fit)` | EAP 信度：var(EAP) / (var(EAP) + mean PSD²) |
| `fit_indices(fit)` | Deviance, pD, DIC, WAIC, LOOIC |
| `pointwise_loglik(fit)` | 抽樣次數 × 人數的邊際 log-likelihood；另有 `fit$loo`、`fit$waic` |
| `posterior_draws(fit, persons = TRUE)` | 抽樣矩陣；`coda::as.mcmc(fit)` 取得 mcmc.list |
| `plot(fit)` | 交叉負荷及 95% 區間（`plot(fit, "loadings")`）；`plot(fit, "trace")` 畫 R̂ 最大的軌跡圖 |
| `ppc(fit)` | 後驗預測檢查：答對率、總分分配、RT 平均數與 SD、條件相依 |
| `convergence(fit)` | Rhat / ESS 摘要（RTMB：最佳化代碼、最大梯度、AGHQ 各輪） |
| `compare(M1 = a, M2 = b)` | 信度、DIC、WAIC、LOOIC、ELPD 差異（`loo_compare`） |
| `jags_code(fit)` | 輸出 JAGS 語法 |
| `rtmb_code(fit)` | RTMB 引擎：輸出獨立的 RTMB 腳本（逐題寫出的聯合密度、Laplace、AGHQ），`source()` 後得到相同的 −2logL |

`birt_quantile` 物件支援 `print()`（各分位數的係數表）、`coef()`、`estimates()`、`reliability()`、`fit_indices()`、`plot()`（β(p) 曲線）。

## 注意事項

- 共變項照原樣使用（請先標準化），不可有遺漏值；指標可以有遺漏（MAR）。
- 邊際概似的積分方式：
  - 第一個潛在變項用 Gauss–Hermite（常態）或 Gauss–Laguerre × Gauss–Hermite（ALD）。
  - 測量層 ALD 這類被積函數有折點的情況，改用密集網格。
  - 第二個潛在變項在指標為常態時解析積分。
  - 已與二維數值積分核對，差距 < 5 × 10⁻⁴。
  - 測量層 ALD 很慢，可先 `fit_indices = FALSE`，之後再 `add_fit_indices(fit, n_draws = 200)`。
- ALD 是 working likelihood，後驗區間可能偏窄（Yang, Wang & He, 2016）。
- 各分位數分開估計，分位數線可能交叉。

## TIMSS 2019 臺灣驗證（N = 631，14 題）

| 模型 | θ 信度 | speed 信度 | LOOIC |
|---|---|---|---|
| M2b（ssp 交叉負荷） | .825 | .827 | 22809.5（與手寫 JAGS 相同） |
| speed ~ θ + 10 個共變項，常態（MR） | .793 | .831 | 23195.5 |
| 同上，ALD p = .5 | .792 | .838 | 23144.9 |

- **極端分位數：** p = .9 時 `speed ~ ability` 為 −0.203（p = .1 / .5 都約為 −0.05）；RTMB 的 SHASH 分配模型（`scripts/24`）得到 −0.176，方向與大小一致。
- 腳本：`scripts/10_birt_timss.R`、`scripts/11_birt_timss_rerun.R`；結果：`results/birt_timss.rds`。
