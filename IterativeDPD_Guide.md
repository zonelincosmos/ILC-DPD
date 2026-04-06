# ILC-DPD (Iterative Learning Control DPD) 技術手冊

---

## 1. 緒論

### 1.1 PA 非線性問題

Power Amplifier (PA) 是無線發射機中不可或缺的元件，負責將信號放大至足夠的功率準位進行傳輸。然而，PA 本質上是非線性元件，當輸入信號振幅接近飽和區域時，會產生兩種主要失真：

- **AM/AM distortion**：輸出振幅相對於理想線性增益產生壓縮。低振幅區域 PA 近似線性，但在高振幅區域輸出逐漸飽和，無法繼續按比例增長。

- **AM/PM distortion**：輸出信號的相位隨輸入振幅而改變。即使輸入信號的相位不變，PA 輸出的相位也會因振幅變化而旋轉。

這兩種非線性效應會造成：
- EVM (Error Vector Magnitude) 升高
- ACLR (Adjacent Channel Leakage Ratio) 惡化，spectral regrowth 增加

### 1.2 DPD 的基本概念

Digital Pre-Distortion (DPD) 的核心思想是：在信號進入 PA 之前，先施加一個與 PA 非線性「互補」的失真，使得信號經過 PA 後，兩者的非線性效應相互抵消，最終輸出近似於線性放大的結果。

```
x(n) ──> [ DPD ] ──> p(n) ──> [ PA ] ──> y(n) ≈ G · x(n)
```

其中 `x(n)` 為原始輸入信號，`p(n)` 為預失真後的信號，`y(n)` 為 PA 輸出，`G` 為期望的 linear gain。若 DPD 能完美補償 PA 非線性，則輸出 `y(n)` 應與 `G · x(n)` 一致。

### 1.3 Model-Free vs. Model-Based DPD

DPD 演算法大致可分為兩類：

**Model-based**（如 memory polynomial、GMP、B-spline）：

建立一個參數化的 predistortion function，透過 least squares 等方式求解係數。訓練完成後，所得的模型（係數或 LUT）可重複使用於相同 PA 上的任意信號。適用於即時通訊系統，能在 digital front-end 持續運作。

**Model-free（ILC-DPD）**：

不假設任何 PA model structure，而是透過反覆量測 PA 輸出，逐 sample 調整 predistorted waveform。每次 iteration 都根據實際量測結果修正波形。最終結果是針對特定輸入信號的最佳預失真波形，而非通用模型。

本文件詳細描述 ILC-DPD 的理論基礎、實作細節、noise handling、memory effect 補償，以及與 B-spline、GMP 等 model-based 方法的比較。

---

## 2. PA Model

範例程式碼 (`IterativeDPD.m`) 提供兩種 PA simulation mode，可透過 `pa_mode` 參數切換：

- **`'memoryless'`**：Rapp AM/AM compression model 加上 AM/PM phase distortion。輸出僅取決於當前輸入 sample `y(n) = f(x(n))`。AM/AM 曲線為單一確定曲線——相同的輸入振幅永遠對應相同的輸出振幅。

- **`'memory'`**：Wiener model——由 4-tap complex FIR filter 串接 memoryless Rapp PA 組成。輸出取決於當前及先前的輸入 samples `y(n) = f(x(n), x(n-1), x(n-2), x(n-3))`。Memory effect 導致 AM/AM 特性出現 **spread**——相同的輸入振幅因 signal history 不同而產生不同的輸出振幅。

所有 PA model 參數均定義於腳本中，可依需求修改。ILC-DPD 演算法本身與 PA model 無關——它只需要能夠發射信號並量測 PA 輸出的能力。

---

## 3. ILC-DPD 演算法

### 3.1 Iterative Learning Control 架構

ILC-DPD 將 predistortion 問題建構為 iterative optimization problem [1][2]：

> 給定期望的 PA 輸出 `y_d(n)`，在不需要 PA parametric model 的條件下，找到 predistorted input `p(n)` 使得 `PA(p(n)) = y_d(n)`。

文獻中存在兩種等價的 iterative update 公式：

**Additive ILC update** [1]：

```
p_{k+1}(n) = p_k(n) + L · ( y_d(n) - y_k(n) )
```

其中 `y_k(n) = PA(p_k(n))` 為第 `k` 次 iteration 的量測輸出，`L` 為 learning gain。此公式直接將 output error 回饋至輸入——若 PA 輸出低於 target，則增加輸入；若高於 target，則減少輸入。收斂條件為 `0 < L < 2/G_ss`，其中 `G_ss` 為 PA small-signal gain [1]。

**Multiplicative update**（本實作採用）：

```
p_{k+1}(n) = p_k(n) · y_d(n) / y_k(n)
```

此公式透過逐 sample 的 complex gain correction，在一次運算中同時補償 AM/AM 和 AM/PM [3]。其物理意義清晰：

- 若 `|y_k(n)| < |y_d(n)|`（PA compression，輸出不足）：比值 > 1，增加 drive
- 若 `|y_k(n)| > |y_d(n)|`（過度驅動）：比值 < 1，降低 drive
- 若 `∠y_k(n) ≠ ∠y_d(n)`（phase shift）：complex division 自動旋轉相位補償

當 `y_k(n) = y_d(n)` 時，比值 = 1，predistorted signal 不再改變——即已達收斂。

### 3.2 演算法詳細步驟

本實作採用 multiplicative update，並引入 learning rate `μ` 與 I/Q averaging 以提高 robustness。

#### Step 0：初始化

```
p_0(n) = x(n)      （identity mapping，尚未施加 predistortion）
```

#### Step 0.5：Saturation Probe

在正式 training 前，需經驗性地確定 PA 的 saturation 特性。此步驟不需要 PA model——純粹透過實際量測完成：

1. 以最大硬體輸入振幅 `r_cap` 驅動 PA
2. 觀測最大輸出 → `P_sat = max(|PA(x_probe)|) × 0.99`
3. 找到對應的輸入振幅 → `r_drive_max`：使 PA 輸出恰好達到 `P_sat` 的輸入振幅

超過 `r_drive_max` 的驅動將無法增加 PA 輸出（已 saturate），因此 predistorted signal 必須被 clamp 在此範圍內。

#### Iteration k = 1, 2, ..., K：

**(a) I/Q-averaged capture**

將當前 predistorted signal `p_k` 通過 PA 傳輸 `N_avg` 次，每次加入獨立的 noise realization，然後取平均：

```
y_bar_k(n) = (1 / N_avg) · Σ_{j=1}^{N_avg} [ PA(p_k(n)) + w_j(n) ]
```

其中 `w_j(n)` 為 AWGN。平均可將 noise variance 降低 `N_avg` 倍。

**(b) Target 計算**

期望輸出為 linear 到 `P_sat`，超過則 hard-clip：

```
y_d(n) = min( G · |x(n)|, P_sat ) / |x(n)| · x(n)
```

振幅超過 `P_sat / G` 的 samples 無法被 linearize——這是 PA 的物理極限，而非演算法的限制。

**(c) Multiplicative correction**

逐 sample complex ratio correction，引入 learning rate `μ` 控制收斂速度：

```
p_{k+1}(n) = p_k(n) · [ (1 - μ) + μ · y_d(n) / y_bar_k(n) ]
```

安全的 complex division 將 magnitude 與 phase 分離處理，避免 divide-by-zero：

```
y_d / y_bar = (|y_d| / |y_bar|) · exp( j · (∠y_d - ∠y_bar) )
```

MATLAB 實作：`y_d ./ max(|y_bar|, epsilon) .* exp(-j * angle(y_bar))`

**(d) Drive clamping**

防止 predistorted signal 超過 PA 的有效輸入範圍：

```
p_{k+1}(n) ← p_{k+1}(n) · min(1, r_drive_max / |p_{k+1}(n)|)
```

超過 `r_drive_max` 的 drive 無法增加 PA 輸出，只會造成浪費並可能損壞元件。

**(e) MSE evaluation**

以 noise-free 的方式評估 DPD 品質（量測的是 predistortion 效果，而非 noise floor）：

```
MSE_k = (1/N) · Σ_{n=1}^{N} | PA(p_k(n)) - G · x(n) |²
```

### 3.3 Convergence

根據 contraction mapping theory [4]，multiplicative iteration 在映射 `p → p · y_d / PA(p)` 為 contraction 時保證收斂。

實際觀察到的收斂行為：

- **Iteration 1–3**：MSE 快速下降（約 10–20 dB），主要非線性效應被修正
- **Iteration 3–10**：逐步改善 saturation boundary 附近的殘餘 distortion
- **Iteration 10 以後**：改善幅度遞減，受限於 `P_sat` hard-clip 和 noise floor

典型結果：memoryless PA 在 20 次 iteration 內可達 20–25 dB 的 MSE 改善。

### 3.4 Memory Effect 的處理

ILC-DPD 相較於 model-based 方法的關鍵優勢：**memory effect 無需顯式建模**。

量測輸出 `y_bar_k(n)` 已經包含了當前 predistorted waveform `p_k` 造成的所有 memory effect。演算法根據實際 PA response 來調整 `p_k(n)`，自然涵蓋了來自鄰近 samples `p_k(n-1), p_k(n-2), ...` 的貢獻。

修改 `p_k(n)` 會透過 memory effect 擾動 `n+1, n+2, ...` 位置的輸出。這些 cross-sample 影響在後續 iteration 中被修正。因此含 memory effect 的 PA 需要更多 iteration 才能收斂（30 次 vs. memoryless 的 20 次），這是因為 sample 間的 coupling 需要更多輪次來穩定。

Model-based 方法（B-spline、GMP）若要處理 memory effect，必須顯式地加入 memory polynomial terms 和 cross-terms，使得複雜度隨 memory depth 增長。

---

## 4. Noise 考量

### 4.1 Model-free 方法的 Noise 敏感性

ILC-DPD 獨立地調整每個 sample（共 `N` 個）。不像 model-based 方法以 `K << N` 個參數 fit `N` 個 data points（overdetermined system 天然具有 smoothing 效果），model-free 方法中 measurement noise 會直接污染 predistorted waveform。這是 model-free 的根本 trade-off。

### 4.2 I/Q Averaging

對 `N_avg` 次獨立 capture 取平均，可將 noise power 降低 `N_avg` 倍：

```
Var(y_bar) = Var(y_single) / N_avg
SNR_eff    = SNR + 10 · log10(N_avg)
```

以 `N_avg = 8` 為例：effective SNR 改善 9 dB（從 30 dB 提升至 39 dB）。

### 4.3 Learning Rate

參數 `μ < 1` 提供額外的 damping，降低單次修正幅度：

- `μ = 1.0`：每次 iteration 完全修正。收斂最快，但最容易受 noise 影響
- `μ = 0.5`：半修正。收斂較慢，但對 noise 更 robust
- `μ = 0.8`：推薦值，在典型 SNR 條件下取得速度與穩定性的平衡

---

## 5. 與 Model-Based DPD 方法的比較

### 5.1 比較方法

三種 DPD 方法在相同的 memoryless PA、相同 waveform、相同 noise 條件下進行評估。詳見 `DPD_Comparison.m`。

**ILC-DPD（本實作）**：逐 sample multiplicative correction，`μ = 0.8`，`N_avg = 8`。

**B-spline DPD**：Complex gain model，degree-3 cubic B-spline basis functions，11 個均勻 break points，K=13 個 basis。透過 WLS (Weighted Least Squares) 訓練。

**GMP DPD（memoryless polynomial）**：Odd-order polynomial gain model，orders {1, 3, 5, 7, 9}，K=5 個 coefficients。採用相同的 WLS training。

### 5.2 實驗結果

Memoryless PA、20 iterations、SNR = 30 dB、N = 30976 samples。

| Method | Final MSE | Improvement | Parameters |
|--------|-----------|-------------|------------|
| No DPD | 4.44e-03 | — | — |
| ILC-DPD | 3.43e-05 | 21.1 dB | N samples |
| GMP (K=5) | 4.74e-05 | 19.7 dB | 5 complex coeff. |
| B-spline (K=13) | 5.50e-05 | 19.1 dB | 13 complex coeff. |

### 5.3 結果分析

ILC-DPD 達到最低 MSE，原因在於每個 sample 均被個別最佳化，不存在 model approximation error。此外，I/Q averaging 提供了 39 dB 的 effective SNR（30 + 9 dB），優於 model-based 方法所使用的單次 30 dB capture。

GMP 優於 B-spline，因為其 5 個 polynomial coefficients 在 noisy 環境下比 13 個 B-spline coefficients 更為 robust。在含 noise 的 training data 中，B-spline 較多的 degrees of freedom 反而導致些許更高的 residual error。

### 5.4 Trade-off 總結

| 面向 | ILC-DPD | B-spline | GMP |
|------|---------|----------|-----|
| MSE | Best | Good | Better |
| 跨信號 reusable | No | Yes | Yes |
| Memory effect 處理 | Implicit | 需顯式 memory terms | 需顯式 memory terms |
| HW deployment | 需儲存完整 waveform | LUT / coefficients | LUT / coefficients |
| Noise robustness | 依賴 I/Q averaging | Model smoothing | Model smoothing |
| 每次 iteration 計算量 | O(N) | O(K²N) | O(K²N) |

### 5.5 適用場景

**ILC-DPD** 適用於 test & calibration 流程——已知 waveform 被重複傳輸的場景，例如：signal generator 校準、type-approval testing、一次性 PA characterization。多家商用 signal generator 廠商已實作此方法。

**B-spline / GMP** 適用於 real-time production transmitter（WiFi、cellular），DPD 必須補償任意 data-bearing signal。訓練後的 model 以 LUT 或 polynomial evaluator 的形式部署於 digital front-end。

---

## 6. 參數表

| 參數 | Symbol | 值 | 說明 |
|------|--------|-----|------|
| Iteration 次數 | K | 20 / 30 | Memoryless / memory PA |
| Learning rate | μ | 0.8 | Correction damping factor |
| I/Q averaging 次數 | N_avg | 8 | 每次 iteration 的 capture 數 |
| Noise SNR | — | 30 dB | ADC measurement noise |
| Target gain | G | 1.0 | 期望 linear gain |
| Probe drive | r_cap | 1.8 | Saturation probe amplitude |
| PA mode | — | `'memoryless'` / `'memory'` | Model 選擇 |

---

## 7. 檔案清單

| 檔案 | 說明 |
|------|------|
| `IterativeDPD.m` | ILC-DPD 實作，支援 memoryless / memory PA 切換 |
| `DPD_Comparison.m` | 三種方法比較：ILC-DPD vs. B-spline vs. GMP |
| `IterativeDPD_Guide.md` | 本文件 |

---

## 8. References

[1] J. Chani-Cahuana, P. N. Landin, C. Fager, and T. Eriksson, "Iterative Learning Control for RF Power Amplifier Linearization," *IEEE Trans. Microw. Theory Techn.*, vol. 64, no. 9, pp. 2778–2789, Sept. 2016.

[2] M. Schoukens, J. Hammenecker, and A. Cooman, "Obtaining the Preinverse of a Power Amplifier Using Iterative Learning Control," *IEEE Trans. Microw. Theory Techn.*, vol. 65, no. 11, pp. 4266–4273, Nov. 2017.

[3] C. Tarver, A. Balatsoukas-Stimming, and J. R. Cavallaro, "Design and Implementation of an Iterative Learning Control Digital Predistortion Algorithm," *Proc. IEEE Asilomar Conf. Signals, Systems, and Computers*, 2018.

[4] US Patent 8,971,829 B2, "Convergence estimation for iterative predistortion factor determination," 2015.
