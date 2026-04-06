# ILC-DPD

Iterative Learning Control Digital Predistortion (ILC-DPD) — a model-free, sample-by-sample waveform predistortion algorithm for power amplifier linearization.

## Overview

ILC-DPD iteratively adjusts each sample of the predistorted waveform based on measured PA output, without requiring a parametric PA model. It naturally handles memory effects through iterative convergence.

See **[IterativeDPD_Guide.md](IterativeDPD_Guide.md)** for the full technical reference including algorithm derivation, noise handling, memory effect analysis, and comparison with model-based approaches.

## Files

| File | Description |
|------|-------------|
| `IterativeDPD.m` | ILC-DPD implementation with memoryless/memory PA selection |
| `DPD_Comparison.m` | Three-method comparison: ILC-DPD vs B-spline vs GMP |
| `IterativeDPD_Guide.md` | Technical reference document |
| `RefX.mat` | Complex baseband test waveform |

## Quick Start

```matlab
cd('path/to/ILC-DPD')

% Run ILC-DPD (edit pa_mode in script: 'memoryless' or 'memory')
IterativeDPD

% Run three-method comparison
DPD_Comparison
```

## Results (SNR = 30 dB)

**Three-method comparison (memoryless PA, 20 iterations):**

| Method | MSE | Improvement |
|--------|-----|-------------|
| ILC-DPD | 3.43e-05 | 21.1 dB |
| GMP (K=5) | 4.74e-05 | 19.7 dB |
| B-spline (K=13) | 5.50e-05 | 19.1 dB |

**Memory PA (Wiener model, 30 iterations):**

| Method | MSE | Improvement |
|--------|-----|-------------|
| ILC-DPD | 3.37e-05 | 26.2 dB |

## References

1. J. Chani-Cahuana et al., "Iterative Learning Control for RF Power Amplifier Linearization," *IEEE Trans. Microw. Theory Techn.*, vol. 64, no. 9, 2016.
2. M. Schoukens et al., "Obtaining the Preinverse of a Power Amplifier Using Iterative Learning Control," *IEEE Trans. Microw. Theory Techn.*, vol. 65, no. 11, 2017.
