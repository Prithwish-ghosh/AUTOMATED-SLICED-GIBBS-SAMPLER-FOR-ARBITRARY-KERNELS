<!DOCTYPE html>

<h1>Automated Sliced Gibbs Sampler for Arbitrary Kernels</h1>

<p>
The arXiv version of the article can be found here.
</p>

<blockquote>
<b>Automated Sliced Gibbs Sampler for Arbitrary Kernels</b><br>
<b>Prithwish Ghosh &amp; Sujit K Ghosh</b><br>
Department of Statistics, North Carolina State University, Raleigh, NC, USA
</blockquote>

<hr>

<h2> Abstract</h2>
We introduce an Automated Sliced Gibbs (ASG) sampling framework for efficient, fully automated Markov chain Monte Carlo (MCMC) inference from arbitrary finite-dimensional probability kernels, including highly non-standard, unnormalized, and non-smooth target densities. The proposed method combines a novel Cauchy-transformation–based effective support estimator with a sequence of slice-driven Gibbs updates, thereby eliminating the need for user-specified truncation bounds, proposal scales, or conditional optimization steps. ASG preserves the invariance and ergodicity properties of classical Gibbs samplers while adapting automatically to multimodal, heavy-tailed, irregular, or otherwise non-Gaussian geometries. A key methodological innovation is its ability to navigate and jump across multiple modes with remarkable ease, producing sequences of nearly independent draws and achieving exceptionally high effective sample size relative to the total post–burn-in sample size.
The ASG framework unifies and extends several strands of adaptive and slice-based MCMC by autonomously detecting high-probability regions, performing localized normalization, and generating updates through bounded univariate slice operations within each Gibbs cycle. We establish the stationarity of the induced Markov chain and show that the automated support calibration does not compromise its convergence properties. Comprehensive numerical experiments demonstrate the broad applicability and robustness of ASG across challenging classes of distributions, including univariate Beta mixtures, multivariate Gaussian mixtures, Rosenbrock and Ackley benchmarks, and non-smooth loss-derived kernels such as LASSO objectives. In all settings, ASG delivers substantially higher effective sample size per second, faster decorrelation, and competitive computational cost compared to standard algorithms such as Random Walk Metropolis–Hastings. By combining automation, geometric adaptability, and theoretical rigor, the ASG sampler provides a scalable, general-purpose solution for sampling from complex, multimodal, and non-standard probability kernels where traditional MCMC methods struggle.

<h2>📌 Overview</h2>

<p>
This repository contains the full implementation and experimental evaluation of the
<b>Automated Sliced Gibbs (ASG)</b> framework for sampling from arbitrary unnormalized kernels.
</p>

<ul>
<li><b>Adaptive effective support estimation</b></li>
<li><b>Slice-based conditional updates</b></li>
<li><b>Automated Gibbs composition</b></li>
<li><b>Diagnostic and time-aware efficiency analysis (ESS and ESS/s)</b></li>
</ul>

<p>
Compared to Metropolis–Hastings (MH), ASG achieves higher statistical efficiency (ESS)
and superior time-aware efficiency (ESS/s).
</p>

<hr>

<h2>✨ Key Features</h2>

<ul>
<li>Works with unnormalized kernels</li>
<li>Automatic dimension inference</li>
<li>No manual truncation bounds</li>
<li>Built-in diagnostics: density, trace, ACF, running mean</li>
<li>Time-based ESS/s evaluation</li>
<li>Tested on:
    <ul>
        <li>Beta mixtures</li>
        <li>Normal mixtures</li>
        <li>Rosenbrock distribution</li>
        <li>Ackley function</li>
    </ul>
</li>
</ul>

<hr>

<h2>⚙️ Execution Order</h2>

<div class="section-box">

<h3>Part I: ASG Sampler Workflow</h3>

<div class="step">
<strong>Step 1:</strong> Load effective support module
<p><code>effective_support_uni_s.R</code></p>
</div>

<div class="step">
<strong>Step 2:</strong> Run core ASG sampler
<p><code>asg_sampler.R</code></p>
</div>

<div class="step">
<strong>Step 3:</strong> Run example scripts
<p>
Examples demonstrate kernel definition, sampling, ESS computation,
and diagnostic plotting.
</p>
</div>

</div>

<div class="section-box">

<h3>Part II: Metropolis–Hastings (MH) Workflow</h3>

<div class="step">
<strong>Step 1:</strong> Initialize MH components
<p><code>mh_initial.R</code></p>
</div>

<div class="step">
<strong>Step 2:</strong> Run MH algorithm
<p><code>mh_algo-2.R</code></p>
</div>

<p>
This provides the baseline MH implementation for comparison
with ASG in terms of ESS and ESS/s.
</p>

</div>

<div class="section-box">

<h3>Part III: Time-Based Efficiency Analysis</h3>

<div class="step">
<strong>Step 1:</strong> Load time-fixed sampling function
<p><code>time fixed function.R</code></p>
</div>

<div class="step">
<strong>Step 2:</strong> Run time-based experiments and generate plots
<p><code>time based analysis examples and plots.R</code></p>
</div>

<p>
These scripts evaluate:
</p>
<ul>
<li>ESS per second (ESS/s)</li>
<li>Time-aware efficiency comparison</li>
<li>Sample-aware efficiency comparison</li>
<li>Performance under fixed run times</li>
</ul>

</div>

<hr>

<h2>📊 Example Results</h2>

<h3>Density Recovery (Univariate)</h3>
<figure>
<p align="center">
<img src="figures/comparisonnnnn.png" width="100%">
</p>
<figcaption align="center">
<b>Figure 1:</b> Density recovery using ASG with effective support estimation.
</figcaption>
</figure>

<h3>Rosenbrock Conditional Density</h3>
<figure>
<p align="center">
<img src="figures/density_overlay.png" width="100%">
</p>
<figcaption align="center">
<b>Figure 2:</b> Conditional density overlay for Rosenbrock kernel.
</figcaption>
</figure>

<h3>Ackley Conditional Density</h3>
<figure>
<p align="center">
<img src="figures/ackley_density.png" width="100%">
</p>
<figcaption align="center">
<b>Figure 3:</b> Conditional density overlay for Ackley kernel.
</figcaption>
</figure>

<hr>

<h3>Time-Based Efficiency (ESS/s)</h3>
<figure>
<p align="center">
<img src="figures/Time_fixed_sample_comparison.png" width="100%">
</p>
<figcaption align="center">
<b>Figure 4:</b> Time-aware efficiency comparison (ESS/s).
</figcaption>
</figure>

<h3>Sample-Based Efficiency (ESS/s)</h3>
<figure>
<p align="center">
<img src="figures/sample_fixed_time_varied_plot_comparison.png" width="100%">
</p>
<figcaption align="center">
<b>Figure 5:</b> Sample-aware efficiency comparison (ESS/s).
</figcaption>
</figure>

<hr>

<p>
For issues, improvements, or collaboration inquiries, please open an issue in this repository.
</p>

</body>
</html>
