# MurineMatchTx

## 1. Introduction

**MurineMatchTx** is an R package for standardized exploratory comparison of human tumors and mouse cancer models at the transcriptomic level. Users define the research question, select the corresponding module and endpoint, and interpret results within the specified input data, ortholog mapping, and analytical parameters.

### Technical contribution

MurineMatchTx connects explicit sample metadata, standardized DEG inputs, audited ortholog mapping, shared threshold objects, input validation, standardized outputs, package-native visualizations, and an executable workflow. It integrates established methods for DEG comparison, uncorrected PCA, TROM scoring, GO enrichment, and exploratory immune and metabolic profiling through documented cross-module data handoffs.

Each module retains its own endpoint. DEG analysis evaluates tumor-versus-normal regulatory concordance, PCA diagnoses global structure and source effects, TROM measures sample-associated gene-set overlap, GO summarizes functional themes, and TME and metabolic modules describe signals derived from human reference signatures. Conflicting module results remain visible and are interpreted against the prespecified research question rather than combined into a universal score.

### Scope of TME and metabolism results

With the current human-reference signatures, TME and metabolism results provide exploratory, method-specific descriptions and remain separate from model-selection endpoints. Cross-species validation requires mouse-specific signatures or orthogonal immune or metabolic evidence.

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
