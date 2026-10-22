# MurineMatchTx

## 1. Introduction

**MurineMatchTx** is an R package designed to address a central challenge in translational cancer research: quantitatively evaluating the **transcriptional similarity between mouse models and human cancer samples**. It provides a streamlined pipeline for cross-species integration of RNA-seq data, enabling researchers to assess the fidelity of preclinical models in mimicking human tumor biology.

By integrating DEG analysis, PCA, TROM scoring, immune profiling, and GO enrichment in a unified workflow, MurineMatchTx helps researchers systematically compare mouse models with human tumor subtypes at the transcriptomic level.

### Scope of TME and metabolism results

When the current human-reference immune-deconvolution or metabolism signatures are applied to mouse samples, their results are exploratory, method-specific descriptions. They do not enter cross-species similarity scoring or model-selection decisions. These results do not establish cross-species validation without mouse-specific signatures or orthogonal immune or metabolic evidence.

### MurineMatchTx package workflow

![workflow](./figure/workflow.png)

---

## 2. Installation

You can install the package directly from GitHub:

```r
# Install devtools if necessary
if (!require("devtools")) install.packages("devtools")

# Install the package
devtools::install_github("wangdian-PKU/MurineMatchTx")
```

---

## 3. User Manual

Detailed guidance of how to use MurineMatchTx could be found in the https://wangdian-pku.github.io/MurineMatchTxManual/.

### CI and analysis environments

The Windows/Linux CI installs dependencies, runs testthat and the executable
vignette, builds the source package, and checks that exact tarball with
`R CMD check --as-cran`. Each run saves commands, logs, `sessionInfo()`, system
information, installed package versions and check results as artifacts. TROM
1.2 is selected from the CRAN archive because it is an archived dependency.
The matrix follows the current R release; it is not a frozen analysis environment.
macOS is deferred until the full dependency chain has a verified installation.

The formal R540_2 analysis project's `renv.lock` belongs to that analysis project.
Creating or updating it requires the stage 5 remote-environment authorization;
this package repository does not create that lockfile or modify the server.

---

## 4. Contact us

We are always open to your questions and encourage discussions here.

You can report potential bugs to the [issues](https://github.com/wangdian-PKU/MurineMatchTx/issues).

Any question can directly contact the author by email: wangdian@bjmu.edu.cn
