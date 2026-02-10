# Automated Sliced Gibbs Sampler for Arbitrary Kernels

This repository contains the code and experiments accompanying the paper:

> **Automated Sliced Gibbs Sampler for Arbitrary Kernels**  
> Prithwish Ghosh, 2025

## 📌 Overview

Statistical inference often requires sampling from complex and intractable probability distributions.  
We introduce an **automated Gibbs sampling framework** that integrates:

- **Adaptive effective support estimation**  
  Automatically identifies and normalizes the support of arbitrary univariate and multivariate kernels.  

- **Slice sampling updates**  
  Uses rejection-based slice sampling to handle nonstandard, multimodal, or highly correlated densities.  

- **Automated Gibbs composition**  
  Extends to multivariate settings without requiring manual bound specification.  

- **Diagnostic and time-based analysis**  
  Includes density overlays, trace plots, running means, autocorrelations, effective sample size (ESS), and ESS per second (ESS/s) for rigorous performance evaluation.  

Compared to Metropolis–Hastings (MH), the proposed method achieves **higher statistical efficiency** (greater ESS) and demonstrates **better time-aware performance** (higher ESS/s under fixed run times).  

---

## ✨ Key Features

- Works directly with **unnormalized kernels**.
- **Automatically infers dimension** by probing kernel evaluations.
- Provides **univariate slice updates** for each conditional distribution.
- Built-in diagnostics:
  - Density plots
  - Trace plots
  - Running mean plots
  - Autocorrelation plots
  - Effective sample size (ESS)
- **Time-based analysis**: evaluates efficiency using ESS per second (ESS/s).
- Tested on:
  - Univariate beta mixtures
  - Multivariate normal mixtures
  - Rosenbrock (“banana-shaped”) distribution
  - Ackley function

---

## 📊 Example Results

### Rosenbrock (“banana-shaped”) distribution

- **MH sampler**: poor mixing, low ESS.  
- **Automated Sliced Gibbs**: adapts to curved geometry, achieves much higher ESS.

<figure>
  <p align="center">
    <img src="figures/rosenbrock_trace.png" alt="Trace plot on Rosenbrock" width="85%">
  </p>
  <figcaption>
    <p align="center">
      <b>Figure 1:</b> Trace plot on the Rosenbrock target. ASG exhibits faster traversal and reduced sticking compared to baseline MH.
    </p>
  </figcaption>
</figure>

<figure>
  <p align="center">
    <img src="figures/rosenbrock_acf.png" alt="ACF on Rosenbrock" width="85%">
  </p>
  <figcaption>
    <p align="center">
      <b>Figure 2:</b> Autocorrelation (ACF) for a representative coordinate. ASG yields lower autocorrelation across lags, indicating better mixing.
    </p>
  </figcaption>
</figure>

---

### Density recovery (univariate conditional example)

<figure>
  <p align="center">
    <img src="figures/univariate_density_overlay.png" alt="Histogram with target overlay" width="85%">
  </p>
  <figcaption>
    <p align="center">
      <b>Figure 3:</b> Histogram of ASG samples with an (initial-conditional) normalized density overlay computed from effective support estimation.
    </p>
  </figcaption>
</figure>

---

### Time-based efficiency (ESS/s)

The proposed sampler maintains higher ESS/s across varying run times, while MH degrades quickly.

<figure>
  <p align="center">
    <img src="figures/ess_per_second.png" alt="ESS per second comparison" width="85%">
  </p>
  <figcaption>
    <p align="center">
      <b>Figure 4:</b> Time-aware efficiency comparison using ESS/s. Higher is better.
    </p>
  </figcaption>
</figure>

---

## 🗂️ Suggested figure folder layout

Place images here:

