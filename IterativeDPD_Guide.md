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

$$x(n) \;\longrightarrow\; \text{DPD} \;\longrightarrow\; p(n) \;\longrightarrow\; \text{PA} \;\longrightarrow\; y(n) \approx G \cdot x(n)$$

where $G$ is the desired linear gain.

### 1.3 Model-Free vs. Model-Based DPD

**Model-based** approaches (memory polynomial, GMP, B-spline) fit a parametric predistortion function $G(|x|)$ and deploy it as a lookup table or polynomial evaluator. The model is reusable for arbitrary input signals on the same PA.

**Model-free (ILC-DPD)** iteratively adjusts the predistorted waveform sample-by-sample, using only measured PA output. No model structure is assumed. The result is a predistorted waveform valid for one specific input signal.

This document describes the ILC-DPD algorithm, its theoretical foundation, noise handling, memory effect compensation, and comparative evaluation against B-spline and GMP approaches.

---

## 2. PA Model

### 2.1 Memoryless PA

The memoryless PA model used throughout this work is the Rapp model with a two-pole AM/PM characteristic.

**AM/AM (Rapp compression model):**

$$\text{amam}(r) = \frac{G_{\text{pa}} \cdot r}{\left(1 + \left(\frac{r}{A_{\text{sat}}}\right)^{2p}\right)^{1/(2p)}}$$

where $r = |x(n)|$ is the instantaneous envelope, $A_{\text{sat}}$ is the saturation amplitude, and $p$ is the smoothness parameter controlling the compression knee.

**AM/PM (two-pole rational model):**

$$\text{ampm}(r) = \left(\frac{\varphi_1 \, r^2}{r_1^2 + r^2} + \frac{\varphi_2 \, r^6}{r_2^6 + r^6}\right) \cdot \frac{\pi}{180}$$

where $\varphi_1, \varphi_2$ are asymptotic phase shifts in degrees, and $r_1, r_2$ are the half-power amplitude points.

**Combined PA response:**

$$y(n) = \frac{\text{amam}(|x(n)|)}{|x(n)|} \cdot x(n) \cdot e^{j \cdot \text{ampm}(|x(n)|)}$$

**Parameters:**

| Symbol | Value | Description |
|--------|-------|-------------|
| $G_{\text{pa}}$ | 1.0 | Small-signal gain |
| $A_{\text{sat}}$ | 0.85 | Saturation amplitude |
| $p$ | 2.5 | Rapp smoothness parameter |
| $\varphi_1$ | 22 deg | Primary AM/PM coefficient |
| $\varphi_2$ | 18 deg | Secondary AM/PM coefficient |
| $r_1$ | 0.70 | Primary AM/PM knee |
| $r_2$ | 0.82 | Secondary AM/PM knee |

### 2.2 Memory PA (Wiener Model)

Real PAs exhibit memory effects: the output $y(n)$ depends on $x(n), x(n\!-\!1), \ldots, x(n\!-\!M)$. Sources include bandwidth-limited matching networks, thermal dynamics, and bias circuit time constants.

The **Wiener model** represents memory as a linear FIR filter preceding the memoryless nonlinearity:

$$x(n) \;\xrightarrow{\text{FIR}}\; \tilde{x}(n) = \sum_{m=0}^{M} h_m \cdot x(n\!-\!m) \;\xrightarrow{\text{NL}}\; y(n) = f_{\text{PA}}(\tilde{x}(n))$$

**FIR coefficients (4-tap, complex):**

| Tap | Value | Magnitude | Physical origin |
|-----|-------|-----------|-----------------|
| $h_0$ | $1.0$ | 100% | Main signal path |
| $h_1$ | $0.15 \, e^{j0.5}$ | 15% | Electrical memory |
| $h_2$ | $-0.10 \, e^{-j0.3}$ | 10% | Impedance mismatch |
| $h_3$ | $0.05 \, e^{j0.7}$ | 5% | Thermal tail |

Memory effects cause the AM/AM characteristic to **spread**: the same input amplitude $|x|$ maps to different output amplitudes depending on signal history.

---

## 3. ILC-DPD Algorithm

### 3.1 Iterative Learning Control Framework

The ILC-DPD formulation treats predistortion as an iterative optimization problem [1][2]:

> Given a desired PA output $y_d(n)$, find a predistorted input $p(n)$ such that $\text{PA}(p(n)) = y_d(n)$, without requiring a parametric model of the PA.

The standard **additive ILC** update [1] is:

$$p_{k+1}(n) = p_k(n) + L \cdot \bigl(y_d(n) - y_k(n)\bigr)$$

where $y_k(n) = \text{PA}(p_k(n))$ is the measured output at iteration $k$, and $L$ is the learning gain.

An equivalent **multiplicative** formulation is:

$$p_{k+1}(n) = p_k(n) \cdot \frac{y_d(n)}{y_k(n)}$$

which applies per-sample complex gain correction. This form naturally handles both AM/AM and AM/PM compensation in a single operation [3].

Both forms converge under appropriate conditions. The multiplicative form requires $|y_k(n)| > 0$ (addressed by numerical safeguards), while the additive form requires $0 < L < 2/G_{\text{ss}}$ for stability [1].

### 3.2 Algorithm Detail

The implementation uses the multiplicative formulation with a learning rate $\mu$ and I/Q averaging.

**Initialization:**

$$p_0(n) = x(n) \quad \text{(identity — no predistortion)}$$

**Saturation probe:** Before training, empirically determine:
- $P_{\text{sat}}$: maximum achievable PA output amplitude (drive PA at $r_{\text{cap}}$, observe output)
- $r_{\text{drive,max}}$: PA input amplitude that produces $P_{\text{sat}}$ (beyond this, PA output cannot increase)

**Iteration $k = 1, 2, \ldots, K$:**

**(a) I/Q-averaged capture.** Transmit $p_k$ through the PA $N_{\text{avg}}$ times with independent noise realizations, and average:

$$\bar{y}_k(n) = \frac{1}{N_{\text{avg}}} \sum_{j=1}^{N_{\text{avg}}} \bigl[\text{PA}(p_k(n)) + w_j(n)\bigr]$$

where $w_j(n)$ is AWGN with variance $\sigma_w^2 = \frac{P_{\text{sig}}}{2 \cdot 10^{\text{SNR}/10}}$ per I/Q component. Averaging reduces noise variance by factor $N_{\text{avg}}$.

**(b) Target computation.** The desired output is linear up to $P_{\text{sat}}$, hard-clipped beyond:

$$y_d(n) = \frac{\min\bigl(G \cdot |x(n)|, \, P_{\text{sat}}\bigr)}{|x(n)|} \cdot x(n)$$

Samples where $G \cdot |x(n)| > P_{\text{sat}}$ cannot be linearized (PA physical limit).

**(c) Multiplicative correction.** Per-sample complex ratio with learning rate $\mu$:

$$p_{k+1}(n) = p_k(n) \cdot \Bigl[(1 - \mu) + \mu \cdot \frac{y_d(n)}{\bar{y}_k(n)}\Bigr]$$

The safe complex division separates magnitude and phase:

$$\frac{y_d}{\bar{y}_k} = \frac{|y_d|}{|\bar{y}_k|} \cdot e^{j(\angle y_d - \angle \bar{y}_k)}$$

implemented as: `y_d ./ max(|y_bar|, epsilon) .* exp(-j * angle(y_bar))`

**(d) Drive clamping.**

$$p_{k+1}(n) \leftarrow p_{k+1}(n) \cdot \min\!\left(1, \; \frac{r_{\text{drive,max}}}{|p_{k+1}(n)|}\right)$$

**(e) MSE evaluation** (noise-free, measures DPD quality):

$$\text{MSE}_k = \frac{1}{N} \sum_{n=1}^{N} \bigl|\text{PA}(p_k(n)) - G \cdot x(n)\bigr|^2$$

### 3.3 Convergence Properties

From contraction mapping theory [4], the multiplicative iteration converges when:

$$\left|\frac{\partial}{\partial p}\left(p \cdot \frac{y_d}{\text{PA}(p)}\right)\right| < 1$$

In practice:
- **Iterations 1--3**: Rapid MSE reduction (~10--20 dB), dominant nonlinearity corrected.
- **Iterations 3--10**: Gradual refinement of residual distortion near the saturation boundary.
- **Beyond 10**: Diminishing returns, limited by $P_{\text{sat}}$ hard-clip and noise floor.

Typical convergence: 20--25 dB MSE improvement within 20 iterations for a memoryless PA.

### 3.4 Memory Effect Handling

A key advantage of ILC-DPD over model-based methods: **memory effects require no explicit modeling**.

The measured output $\bar{y}_k(n)$ already incorporates all memory effects from the current predistorted waveform $p_k$. The per-sample correction adjusts $p_k(n)$ based on the actual PA response, including contributions from neighboring samples $p_k(n\!-\!1), p_k(n\!-\!2), \ldots$

Modifying $p_k(n)$ perturbs the output at $n\!+\!1, n\!+\!2, \ldots$ through the memory taps. These cross-sample effects are corrected in subsequent iterations. More iterations are needed for convergence (30 vs. 20 for memoryless PA) due to this coupling.

Model-based methods (B-spline, GMP) must explicitly add memory polynomial terms $x(n\!-\!q) \cdot |x(n\!-\!q)|^{p-1}$ and cross-terms, increasing complexity exponentially with memory depth.

---

## 4. Noise Considerations

### 4.1 Noise Sensitivity

ILC-DPD adjusts each of $N$ samples independently. Without smoothing from a model fit ($K \ll N$ parameters), measurement noise directly corrupts the predistorted waveform. This is the fundamental trade-off of model-free operation.

### 4.2 I/Q Averaging

Averaging $N_{\text{avg}}$ independent captures reduces noise power by factor $N_{\text{avg}}$:

$$\text{Var}(\bar{y}_k) = \frac{\text{Var}(y_{\text{single}})}{N_{\text{avg}}}, \qquad \text{SNR}_{\text{eff}} = \text{SNR} + 10 \log_{10}(N_{\text{avg}})$$

With $N_{\text{avg}} = 8$: effective SNR improves by 9 dB.

### 4.3 Learning Rate

The parameter $\mu < 1$ provides additional damping:

$$p_{k+1} = p_k \cdot \bigl[(1 - \mu) + \mu \cdot y_d / \bar{y}_k\bigr]$$

At $\mu = 1$: full correction per iteration (fastest, most noise-sensitive).
At $\mu = 0.5$: half correction (slower convergence, more robust to noise).
Recommended: $\mu = 0.8$ for typical SNR conditions.

---

## 5. Comparison with Model-Based DPD

### 5.1 Methods Compared

Three DPD approaches are evaluated on the same memoryless PA, waveform, and noise conditions.

**ILC-DPD (this work):**
Per-sample multiplicative correction, $\mu = 0.8$, $N_{\text{avg}} = 8$.

**B-spline DPD:**
Complex gain model $G(r) = \sum_{k=1}^{K} c_k \, B_k(r)$ with degree-3 cubic B-splines, 11 uniform break points on $[0, 1]$, $K = 13$ basis functions. Trained via weighted least squares (WLS) with amplitude-weighted cost $w(r) = 1 + 20 \cdot r^4$ and Tikhonov regularization $\lambda = 10^{-7}$.

**GMP DPD (memoryless polynomial):**
Complex gain model $G(r) = \sum_{k} c_k \, r^{p_k - 1}$ with odd orders $p_k \in \{1, 3, 5, 7, 9\}$, $K = 5$ coefficients. Same WLS training.

### 5.2 Experimental Results

Memoryless Rapp PA, 20 iterations, SNR = 30 dB, $N = 30976$ samples.

| Method | Final MSE | Improvement | Parameters |
|--------|-----------|-------------|------------|
| No DPD | 4.44e-03 | --- | --- |
| ILC-DPD | 3.43e-05 | 21.1 dB | $N$ samples |
| GMP ($K\!=\!5$) | 4.74e-05 | 19.7 dB | 5 complex coeff. |
| B-spline ($K\!=\!13$) | 5.50e-05 | 19.1 dB | 13 complex coeff. |

### 5.3 Analysis

ILC-DPD achieves the lowest MSE because each sample is individually optimized with no model approximation error. I/Q averaging provides an effective SNR of 39 dB ($30 + 10\log_{10}(8) \approx 39$ dB), exceeding the single-capture SNR of 30 dB available to the model-based methods.

GMP outperforms B-spline because its 5 polynomial coefficients are more robust to noise than 13 B-spline coefficients. With noisy training data, the additional degrees of freedom in the B-spline model can lead to slightly higher residual error.

### 5.4 Trade-off Summary

| Aspect | ILC-DPD | B-spline | GMP |
|--------|---------|----------|-----|
| MSE | Best | Good | Better |
| Reusable across signals | No | Yes | Yes |
| Memory effect handling | Implicit | Requires memory taps | Requires memory taps |
| HW deployment | Full waveform storage | LUT / coefficients | LUT / coefficients |
| Noise robustness | I/Q averaging | Model smoothing | Model smoothing |
| Computation per iter | $O(N)$ | $O(K^2 N)$ | $O(K^2 N)$ |

### 5.5 Application Scenarios

**ILC-DPD** is suited for test and calibration workflows where a known waveform is transmitted repeatedly: signal generator calibration, type-approval testing, one-time PA characterization. Commercial signal generators from multiple vendors implement this approach.

**B-spline / GMP** are suited for real-time production transmitters (WiFi, cellular) where the DPD must compensate arbitrary data-bearing signals. The trained model is deployed as a lookup table or polynomial evaluator in the digital front-end.

---

## 6. Parameters

| Parameter | Symbol | Value | Description |
|-----------|--------|-------|-------------|
| Iterations | $K$ | 20 / 30 | Memoryless / memory PA |
| Learning rate | $\mu$ | 0.8 | Correction damping factor |
| I/Q averages | $N_{\text{avg}}$ | 8 | Captures per iteration |
| Noise SNR | --- | 30 dB | ADC measurement noise |
| Target gain | $G$ | 1.0 | Desired linear gain |
| Probe drive | $r_{\text{cap}}$ | 1.8 | Saturation probe amplitude |
| PA mode | --- | `'memoryless'` / `'memory'` | Model selection |

---

## 7. File Inventory

| File | Description |
|------|-------------|
| `IterativeDPD.m` | ILC-DPD implementation with memoryless/memory PA selection |
| `DPD_Comparison.m` | Three-method comparison: ILC-DPD vs. B-spline vs. GMP |
| `IterativeDPD_Guide.md` | This document |

---

## 8. References

[1] J. Chani-Cahuana, P. N. Landin, C. Fager, and T. Eriksson, "Iterative Learning Control for RF Power Amplifier Linearization," *IEEE Trans. Microw. Theory Techn.*, vol. 64, no. 9, pp. 2778--2789, Sept. 2016.

[2] M. Schoukens, J. Hammenecker, and A. Cooman, "Obtaining the Preinverse of a Power Amplifier Using Iterative Learning Control," *IEEE Trans. Microw. Theory Techn.*, vol. 65, no. 11, pp. 4266--4273, Nov. 2017.

[3] C. Tarver, A. Balatsoukas-Stimming, and J. R. Cavallaro, "Design and Implementation of an Iterative Learning Control Digital Predistortion Algorithm," *Proc. IEEE Asilomar Conf. Signals, Systems, and Computers*, 2018.

[4] US Patent 8,971,829 B2, "Convergence estimation for iterative predistortion factor determination," 2015.
