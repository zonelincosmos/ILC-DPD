# Iterative Learning Control Digital Predistortion (ILC-DPD)

## Technical Reference

---

## 1. Introduction

### 1.1 Power Amplifier Nonlinearity

A power amplifier (PA) is inherently nonlinear. As the input signal amplitude approaches saturation, the PA exhibits:

- **AM/AM distortion**: output amplitude compresses relative to the ideal linear gain.
- **AM/PM distortion**: output phase rotates as a function of input amplitude.

These effects degrade modulation accuracy (EVM) and generate spectral regrowth (ACLR).

### 1.2 Digital Predistortion

Digital predistortion (DPD) compensates PA nonlinearity by applying an inverse distortion to the input signal before the PA:

```
x(n) ──> [ DPD ] ──> p(n) ──> [ PA ] ──> y(n) ≈ G · x(n)
```

where `G` is the desired linear gain.

### 1.3 Model-Free vs. Model-Based DPD

**Model-based** approaches (memory polynomial, GMP, B-spline) fit a parametric predistortion function and deploy it as a lookup table or polynomial evaluator. The model is reusable for arbitrary input signals on the same PA.

**Model-free (ILC-DPD)** iteratively adjusts the predistorted waveform sample-by-sample, using only measured PA output. No model structure is assumed. The result is a predistorted waveform valid for one specific input signal.

This document describes the ILC-DPD algorithm, its theoretical foundation, noise handling, memory effect compensation, and comparative evaluation against B-spline and GMP approaches.

---

## 2. PA Model

The example code (`IterativeDPD.m`) provides two PA simulation modes, selectable via the `pa_mode` parameter:

- **`'memoryless'`** — Rapp AM/AM compression + AM/PM phase distortion. Output depends only on the current input sample: `y(n) = f(x(n))`.

- **`'memory'`** — Wiener model: a 4-tap complex FIR filter followed by the memoryless Rapp PA. Output depends on the current and previous input samples: `y(n) = f(x(n), x(n-1), x(n-2), x(n-3))`. Memory effects cause the AM/AM characteristic to **spread** — the same input amplitude produces different output amplitudes depending on signal history.

All PA model parameters are defined in the script and can be modified as needed. The ILC-DPD algorithm itself is PA-model-agnostic — it only requires the ability to transmit a signal and measure the PA output.

---

## 3. ILC-DPD Algorithm

### 3.1 Iterative Learning Control Framework

The ILC-DPD formulation treats predistortion as an iterative optimization problem [1][2]:

> Given a desired PA output `y_d(n)`, find a predistorted input `p(n)` such that `PA(p(n)) = y_d(n)`, without requiring a parametric model of the PA.

The standard **additive ILC** update [1] is:

```
p_{k+1}(n) = p_k(n) + L · ( y_d(n) - y_k(n) )
```

where `y_k(n) = PA(p_k(n))` is the measured output at iteration `k`, and `L` is the learning gain.

An equivalent **multiplicative** formulation is:

```
p_{k+1}(n) = p_k(n) · y_d(n) / y_k(n)
```

which applies per-sample complex gain correction. This form naturally handles both AM/AM and AM/PM compensation in a single operation [3].

Both forms converge under appropriate conditions. The multiplicative form requires `|y_k(n)| > 0` (addressed by numerical safeguards), while the additive form requires `0 < L < 2/G_ss` for stability [1].

### 3.2 Algorithm Detail

The implementation uses the multiplicative formulation with a learning rate `μ` and I/Q averaging.

**Initialization:**

```
p_0(n) = x(n)      (identity — no predistortion)
```

**Saturation probe:** Before training, empirically determine:
- `P_sat`: maximum achievable PA output amplitude (drive PA at `r_cap`, observe output)
- `r_drive_max`: PA input amplitude that produces `P_sat` (beyond this, PA output cannot increase)

**Iteration k = 1, 2, ..., K:**

**(a) I/Q-averaged capture.** Transmit `p_k` through the PA `N_avg` times with independent noise realizations, and average:

```
y_bar_k(n) = (1 / N_avg) · Σ_{j=1}^{N_avg} [ PA(p_k(n)) + w_j(n) ]
```

where `w_j(n)` is AWGN. Averaging reduces noise variance by factor `N_avg`.

**(b) Target computation.** The desired output is linear up to `P_sat`, hard-clipped beyond:

```
y_d(n) = min( G · |x(n)|, P_sat ) / |x(n)| · x(n)
```

Samples where `G · |x(n)| > P_sat` cannot be linearized (PA physical limit).

**(c) Multiplicative correction.** Per-sample complex ratio with learning rate `μ`:

```
p_{k+1}(n) = p_k(n) · [ (1 - μ) + μ · y_d(n) / y_bar_k(n) ]
```

The safe complex division separates magnitude and phase:

```
y_d / y_bar = (|y_d| / |y_bar|) · exp( j · (∠y_d − ∠y_bar) )
```

**(d) Drive clamping.**

```
p_{k+1}(n) ← p_{k+1}(n) · min(1, r_drive_max / |p_{k+1}(n)|)
```

**(e) MSE evaluation** (noise-free, measures DPD quality):

```
MSE_k = (1/N) · Σ_{n=1}^{N} | PA(p_k(n)) − G · x(n) |²
```

### 3.3 Convergence Properties

From contraction mapping theory [4], the multiplicative iteration converges when the map `p → p · y_d / PA(p)` is a contraction.

In practice:
- **Iterations 1–3**: Rapid MSE reduction (~10–20 dB), dominant nonlinearity corrected.
- **Iterations 3–10**: Gradual refinement of residual distortion near the saturation boundary.
- **Beyond 10**: Diminishing returns, limited by `P_sat` hard-clip and noise floor.

Typical convergence: 20–25 dB MSE improvement within 20 iterations for a memoryless PA.

### 3.4 Memory Effect Handling

A key advantage of ILC-DPD over model-based methods: **memory effects require no explicit modeling**.

The measured output `y_bar_k(n)` already incorporates all memory effects from the current predistorted waveform `p_k`. The per-sample correction adjusts `p_k(n)` based on the actual PA response, including contributions from neighboring samples `p_k(n-1), p_k(n-2), ...`

Modifying `p_k(n)` perturbs the output at `n+1, n+2, ...` through the memory taps. These cross-sample effects are corrected in subsequent iterations. More iterations are needed for convergence (30 vs. 20 for memoryless PA) due to this coupling.

Model-based methods (B-spline, GMP) must explicitly add memory polynomial terms and cross-terms, increasing complexity with memory depth.

---

## 4. Noise Considerations

### 4.1 Noise Sensitivity

ILC-DPD adjusts each of `N` samples independently. Without smoothing from a model fit (`K << N` parameters), measurement noise directly corrupts the predistorted waveform. This is the fundamental trade-off of model-free operation.

### 4.2 I/Q Averaging

Averaging `N_avg` independent captures reduces noise power by factor `N_avg`:

```
Var(y_bar) = Var(y_single) / N_avg
SNR_eff    = SNR + 10·log10(N_avg)
```

With `N_avg = 8`: effective SNR improves by 9 dB.

### 4.3 Learning Rate

The parameter `μ < 1` provides additional damping:

- `μ = 1.0`: full correction per iteration (fastest, most noise-sensitive).
- `μ = 0.5`: half correction (slower convergence, more robust to noise).
- `μ = 0.8`: recommended balance for typical SNR conditions.

---

## 5. Comparison with Model-Based DPD

### 5.1 Methods Compared

Three DPD approaches are evaluated on the same memoryless PA, waveform, and noise conditions. See `DPD_Comparison.m`.

**ILC-DPD (this work):**
Per-sample multiplicative correction, `μ = 0.8`, `N_avg = 8`.

**B-spline DPD:**
Complex gain model with degree-3 cubic B-splines, 11 uniform break points, K=13 basis functions. Trained via weighted least squares (WLS).

**GMP DPD (memoryless polynomial):**
Odd-order polynomial gain model, orders {1, 3, 5, 7, 9}, K=5 coefficients. Same WLS training.

### 5.2 Experimental Results

Memoryless PA, 20 iterations, SNR = 30 dB, N = 30976 samples.

| Method | Final MSE | Improvement | Parameters |
|--------|-----------|-------------|------------|
| No DPD | 4.44e-03 | — | — |
| ILC-DPD | 3.43e-05 | 21.1 dB | N samples |
| GMP (K=5) | 4.74e-05 | 19.7 dB | 5 complex coeff. |
| B-spline (K=13) | 5.50e-05 | 19.1 dB | 13 complex coeff. |

### 5.3 Analysis

ILC-DPD achieves the lowest MSE because each sample is individually optimized with no model approximation error. I/Q averaging provides an effective SNR of 39 dB (30 + 9 dB from N_avg=8), exceeding the single-capture SNR of 30 dB available to the model-based methods.

GMP outperforms B-spline because its 5 polynomial coefficients are more robust to noise than 13 B-spline coefficients. With noisy training data, the additional degrees of freedom in the B-spline model can lead to slightly higher residual error.

### 5.4 Trade-off Summary

| Aspect | ILC-DPD | B-spline | GMP |
|--------|---------|----------|-----|
| MSE | Best | Good | Better |
| Reusable across signals | No | Yes | Yes |
| Memory effect handling | Implicit | Requires memory taps | Requires memory taps |
| HW deployment | Full waveform storage | LUT / coefficients | LUT / coefficients |
| Noise robustness | I/Q averaging | Model smoothing | Model smoothing |
| Computation per iter | O(N) | O(K²N) | O(K²N) |

### 5.5 Application Scenarios

**ILC-DPD** is suited for test and calibration workflows where a known waveform is transmitted repeatedly: signal generator calibration, type-approval testing, one-time PA characterization. Commercial signal generators from multiple vendors implement this approach.

**B-spline / GMP** are suited for real-time production transmitters (WiFi, cellular) where the DPD must compensate arbitrary data-bearing signals. The trained model is deployed as a lookup table or polynomial evaluator in the digital front-end.

---

## 6. Parameters

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Iterations | K | 20 / 30 | Memoryless / memory PA |
| Learning rate | μ | 0.8 | Correction damping factor |
| I/Q averages | N_avg | 8 | Captures per iteration |
| Noise SNR | — | 30 dB | ADC measurement noise |
| Target gain | G | 1.0 | Desired linear gain |
| Probe drive | r_cap | 1.8 | Saturation probe amplitude |
| PA mode | — | `'memoryless'` / `'memory'` | Model selection |

---

## 7. File Inventory

| File | Description |
|------|-------------|
| `IterativeDPD.m` | ILC-DPD implementation with memoryless/memory PA selection |
| `DPD_Comparison.m` | Three-method comparison: ILC-DPD vs. B-spline vs. GMP |
| `IterativeDPD_Guide.md` | This document |

---

## 8. References

[1] J. Chani-Cahuana, P. N. Landin, C. Fager, and T. Eriksson, "Iterative Learning Control for RF Power Amplifier Linearization," *IEEE Trans. Microw. Theory Techn.*, vol. 64, no. 9, pp. 2778–2789, Sept. 2016.

[2] M. Schoukens, J. Hammenecker, and A. Cooman, "Obtaining the Preinverse of a Power Amplifier Using Iterative Learning Control," *IEEE Trans. Microw. Theory Techn.*, vol. 65, no. 11, pp. 4266–4273, Nov. 2017.

[3] C. Tarver, A. Balatsoukas-Stimming, and J. R. Cavallaro, "Design and Implementation of an Iterative Learning Control Digital Predistortion Algorithm," *Proc. IEEE Asilomar Conf. Signals, Systems, and Computers*, 2018.

[4] US Patent 8,971,829 B2, "Convergence estimation for iterative predistortion factor determination," 2015.
