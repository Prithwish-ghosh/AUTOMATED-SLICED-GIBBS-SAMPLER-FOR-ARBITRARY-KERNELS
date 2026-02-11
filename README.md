

<h1>Automated Sliced Gibbs Sampler for Arbitrary Kernels</h1>

<p>
The arXiv version of the article can be found here.
</p>

<p>
This repository contains the code and experiments accompanying the paper:
</p>

<blockquote>
<b>Automated Sliced Gibbs Sampler for Arbitrary Kernels</b><br>
<b>Prithwish Ghosh &amp; Sujit K Ghosh</b><br>
Department of Statistics, North Carolina State University, Raleigh, NC, USA
</blockquote>

<hr>

<h2>📌 Overview</h2>

<p>
Statistical inference often requires sampling from complex and intractable probability distributions.
We introduce an <b>automated Gibbs sampling framework</b> that integrates:
</p>

<ul>
<li><b>Adaptive effective support estimation</b><br>
Automatically identifies and normalizes the support of arbitrary univariate and multivariate kernels.</li>

<li><b>Slice sampling updates</b><br>
Uses rejection-based slice sampling to handle nonstandard, multimodal, or highly correlated densities.</li>

<li><b>Automated Gibbs composition</b><br>
Extends to multivariate settings without requiring manual bound specification.</li>

<li><b>Diagnostic and time-based analysis</b><br>
Includes density overlays, trace plots, running means, autocorrelations, effective sample size (ESS), and ESS per second (ESS/s).</li>
</ul>

<p>
Compared to Metropolis–Hastings (MH), the proposed method achieves 
<b>higher statistical efficiency</b> (greater ESS) and demonstrates 
<b>better time-aware performance</b> (higher ESS/s under fixed run times).
</p>

<hr>

<h2>✨ Key Features</h2>

<ul>
<li>Works directly with <b>unnormalized kernels</b>.</li>
<li><b>Automatically infers dimension</b> by probing kernel evaluations.</li>
<li>Provides <b>univariate slice updates</b> for each conditional distribution.</li>
<li>Built-in diagnostics:
    <ul>
        <li>Density plots</li>
        <li>Trace plots</li>
        <li>Running mean plots</li>
        <li>Autocorrelation plots</li>
        <li>Effective sample size (ESS)</li>
    </ul>
</li>
<li><b>Time-based analysis</b>: evaluates efficiency using ESS per second (ESS/s).</li>
<li>Tested on:
    <ul>
        <li>Univariate beta mixtures</li>
        <li>Multivariate normal mixtures</li>
        <li>Rosenbrock (“banana-shaped”) distribution</li>
        <li>Ackley function</li>
    </ul>
</li>
</ul>

<hr>

<h2>⚙️ Execution Order</h2>

<div class="step">
<h3>Step 1: Load Effective Support Module</h3>
<p>First run:</p>
<p><code>effective_support_uni_s.R</code></p>
<p>
This script estimates the effective support of the target distribution.
It automatically calibrates the sampling region and removes the need for manual truncation bounds.
</p>
</div>

<div class="step">
<h3>Step 2: Run the ASG Sampler</h3>
<p>Then run:</p>
<p><code>asg_sampler.R</code></p>
<p>
This file contains the core Automated Sliced Gibbs implementation.
It performs slice-driven updates and generates near-independent posterior samples.
</p>
</div>

<div class="step">
<h3>Step 3: Run the Examples</h3>
<p>
After loading both modules, run the example scripts to reproduce experiments and figures.
</p>
<ul>
<li>Define a target kernel</li>
<li>Initialize parameters</li>
<li>Generate posterior samples</li>
<li>Compute ESS and ESS/s</li>
<li>Generate diagnostic plots</li>
</ul>
</div>

<hr>

<h2>📊 Example Results</h2>

<h3>Density Recovery (Univariate Example)</h3>

<figure>
<p align="center">
<img src="figures/comparisonnnnn.png" width="85%">
</p>
<figcaption align="center">
<b>Figure 1:</b> Sample histogram of ASG samples with density overlay computed from effective support estimation.
</figcaption>
</figure>

<h3>Density Recovery (Rosenbrock Kernel)</h3>

<figure>
<p align="center">
<img src="figures/density_overlay.png" width="85%">
</p>
<figcaption align="center">
<b>Figure 2:</b> Sample histogram with density overlay for Rosenbrock conditional.
</figcaption>
</figure>

<h3>Density Recovery (Ackley Kernel)</h3>

<figure>
<p align="center">
<img src="figures/ackley_density.png" width="85%">
</p>
<figcaption align="center">
<b>Figure 3:</b> Sample histogram with density overlay for Ackley conditional.
</figcaption>
</figure>

<hr>

<h3>Time-Based Efficiency (ESS/s)</h3>

<figure>
<p align="center">
<img src="figures/Time_fixed_sample_comparison.png" width="85%">
</p>
<figcaption align="center">
<b>Figure 4:</b> Time-aware efficiency comparison using ESS/s.
</figcaption>
</figure>

<h3>Sample-Based Efficiency (ESS/s)</h3>

<figure>
<p align="center">
<img src="figures/sample_fixed_time_varied_plot_comparison.png" width="85%">
</p>
<figcaption align="center">
<b>Figure 5:</b> Sample-aware efficiency comparison using ESS/s.
</figcaption>
</figure>

<hr>

<p>
For questions, suggestions, or collaborations, please open an issue in this repository.
</p>

</body>
</html>
