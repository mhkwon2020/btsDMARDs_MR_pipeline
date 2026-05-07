# btsDMARDs MR Pipeline (Version 1)

## Overview

This repository provides a fully reproducible Mendelian randomisation (MR) pipeline for evaluating drug target validity in rheumatoid arthritis (RA) using multi-layer genetic evidence, including eQTL, pQTL, colocalisation, and triangulation frameworks.

The pipeline is **protocol-driven, pre-specified, and designed for publication-grade reproducibility**, enabling transparent and auditable causal inference in genetic epidemiology.

---

## Scientific Objective

- Evaluate causal effects of 18 genetically proxied bDMARD/tsDMARD targets on RA risk
- Integrate multi-omics evidence (cis-eQTL, cis-pQTL, colocalisation)
- Identify mechanistic failure modes in MR inference
- Provide a structured triangulation framework for drug target validation

---

## Quick Start

### 1. Prerequisites

- R ≥ 4.0
- PLINK 1.9 ([download](https://www.cog-genomics.org/plink/)) — optional if using pre-clumped instruments (see [LD Clumping](#ld-clumping))
- OpenGWAS JWT token — **required** for outcome data access (see [OpenGWAS API Token](#opengwas-api-token))
- External data files (see [External Data Configuration](#external-data-configuration))

### 2. Set up environment variables

Copy `.Renviron.example` to `.Renviron` in the project root and fill in your local paths:

```bash
cp .Renviron.example .Renviron
```

Then edit `.Renviron`:

```
PLINK_BIN=/path/to/plink
LD_REF_EUR=/path/to/1000G/EUR
LD_REF_EAS=/path/to/1000G/EAS
INPUT_DATA_PATH=/path/to/eqtl/exposure_eqtlgen_targets18_betaSE.txt
PQTL_DATA_DIR=/path/to/pqtl
OPENGWAS_JWT=<your_token>
```

> **Note:** `.Renviron` is listed in `.gitignore` and will never be committed to version control.

### 3. Run the pipeline

Open R, set the working directory to the project root, then:

```r
# Load environment variables (if not auto-loaded at R startup)
readRenviron(".Renviron")

# Run full pipeline
source("main_fullstudy.R")
```

> **Tip (Windows):** `.Renviron` is auto-loaded when R starts from the project directory.  
> If it is not picked up automatically (e.g., when using a system-level `.Renviron`), call  
> `readRenviron("C:/path/to/project/.Renviron")` before sourcing.

---

## OpenGWAS API Token

The pipeline queries the [OpenGWAS platform](https://gwas.mrcieu.ac.uk/) for outcome GWAS summary statistics and for the failure mode diagnostic API fetch layer. Both require a **free JWT authentication token**, obtainable from [https://api.opengwas.io/](https://api.opengwas.io/).

The `OPENGWAS_JWT` environment variable is used by `04_outcome_harmonise.R` (outcome data extraction) and `15_failure_mode_diagnostic.R` (external API annotation for diagnostic categories 2, 4, and 6). The token must be set in `.Renviron` before running the pipeline.

> **For collaborators:** Each user must obtain and configure their own JWT token. Tokens are personal credentials and must **never** be committed to version control or shared in any configuration file.

---

## LD Clumping

Pre-clumped instrument files for all 18 targets are included in `cache/clumped_instruments/`.  
The pipeline checks this cache before calling PLINK, so **PLINK is not required** if you are reproducing the exact analysis with the same protocol parameters (p = 5×10⁻⁸, r² = 0.001, kb = 10,000, EUR ancestry).

If you change these parameters or add new targets, PLINK and the 1000G LD reference panels will be needed.

---

## Pre-extracted ARIES mQTL Data

A pre-extracted ARIES mQTL dataset for all 18 drug targets is included in `cache/aries_mqtl_17genes.csv`.  
This file was extracted from the `MRInstruments` R package (`MRInstruments::aries_mqtl`) and filtered to the 18 target genes. It is used by the Category 6 (post-transcriptional discordance) failure mode diagnostic.

The pipeline (`R/16_failure_mode_cat6.R`) searches for ARIES data in the following order:
1. `<run_dir>/results/aries_mqtl_17genes.csv` — run-specific output (if present)
2. `cache/aries_mqtl_17genes.csv` — pre-extracted file included in this repository
3. Per-gene `.rds` files in `cache/failmode_api/` — generated at runtime if API fetch is enabled

Because the pre-extracted CSV is committed to the repository, **the `MRInstruments` package is not required** for standard reproduction.  
If you need to regenerate the CSV (e.g., for an updated target list), install and run:

```r
remotes::install_github("MRCIEU/MRInstruments")
library(MRInstruments)
data(aries_mqtl)
targets_upper <- toupper(c("tnf","lta","fcgr1a","fcgr2a","fcgr2b","fcgr2c",
                            "fcgr3a","fcgr3b","cd80","cd86","ctla4","ms4a1",
                            "il6r","il1r1","jak1","jak2","jak3","tyk2"))
aries_sub <- aries_mqtl[aries_mqtl$gene %in% targets_upper, ]
write.csv(aries_sub, "cache/aries_mqtl_17genes.csv", row.names = FALSE)
```

---

## Reproducibility Architecture

### Pipeline Flow

```
Protocol (YAML)
      │
      ▼
Input Readers
      │
      ▼
Instrument Selection (eQTL)
      │
      ▼
LD Clumping & Harmonisation
      │
      ▼
MR Estimation
      │
 ┌────┴────┐
 ▼         ▼
Colocalisation   pQTL Validation
      │         │
      └────┬────┘
           ▼
     Classification
           │
           ▼
     Triangulation
           │
           ▼
 Failure Mode Diagnostic
           │
           ▼
        Reporting
```

### Versioning Strategy

| Component         | Versioning                  |
| ----------------- | --------------------------- |
| R scripts         | Stable (no version suffix)  |
| Protocol (YAML)   | Versioned (v1)              |
| Run configuration | Versioned                   |
| Output            | Timestamped                 |
| Analysis pipeline | Version 1 (current release) |

---

## Project Structure

```
project_root/
  README.md
  main_fullstudy.R
  .Renviron.example          # Template — copy to .Renviron and fill in paths
  LICENSE

  R/
    00_utils.R               # Core utilities (logging, env path expansion)
    01_protocol_checks.R     # Protocol validation
    02_input_readers.R       # Data ingestion helpers
    03_instrument_builder.R  # Instrument selection & LD clumping
    04_outcome_harmonise.R   # Outcome extraction & harmonisation
    05_sample_overlap_check.R
    06_pqtl_loader.R
    07_mr_estimation.R
    08_pqtl_validation.R
    09_colocalisation.R
    10_classification.R
    11_reporting_tables.R
    12_reporting_figures.R
    13_triangulation.R
    14_additional_figures.R
    15_failure_mode_diagnostic.R
    16_failure_mode_cat6.R

  configs/
    protocol_v1.yaml              # Pre-specified analysis protocol
    run_full_study_v1.yaml        # Run-level execution settings
    targets.yaml                  # 18 drug target metadata
    pqtl_sources.yaml             # pQTL data source manifest
    sample_overlap_metadata.yaml  # Sample overlap reference
    failure_mode_diagnostic.yaml  # Failure mode classifier config

  cache/
    clumped_instruments/          # Pre-clumped instruments (EUR, p=5e-8)
                                  # Allows running without PLINK
    aries_mqtl_17genes.csv        # Pre-extracted ARIES mQTL for 18 targets
                                  # Allows Category 6 analysis without MRInstruments
    failmode_api/                 # Per-gene API fetch cache (gitignored; runtime-generated)

  runs/
    <timestamp>/                  # Auto-generated per run
      artifacts/
        config_snapshot.yaml
        session_info.txt
```

---

## External Data Configuration

Large external data files are **not included** in this repository.  
Paths are resolved at runtime via environment variables (see `.Renviron.example`).

| Variable | Description |
|---|---|
| `INPUT_DATA_PATH` | eQTLGen cis-eQTL exposure file (`.txt`) |
| `PQTL_DATA_DIR` | Directory containing pQTL files |
| `LD_REF_EUR` | 1000G EUR PLINK bfile prefix (`.bed/.bim/.fam`) |
| `LD_REF_EAS` | 1000G EAS PLINK bfile prefix |
| `PLINK_BIN` | Path to PLINK 1.9 binary |
| `OPENGWAS_JWT` | JWT token for OpenGWAS API (see [OpenGWAS API Token](#opengwas-api-token)) |

### Required pQTL files

```
${PQTL_DATA_DIR}/ukb_ppp_targets18.tsv
${PQTL_DATA_DIR}/decode_main_clumped_mr_ready.tsv.gz
```

### Required eQTL file

```
${INPUT_DATA_PATH}   # e.g., exposure_eqtlgen_targets18_betaSE.txt
```

---

## Methods ↔ Code Mapping

| Methods Component             | Code Module                    |
| ----------------------------- | ------------------------------ |
| Study protocol definition     | `01_protocol_checks.R`         |
| Data ingestion                | `02_input_readers.R`           |
| Instrument selection          | `03_instrument_builder.R`      |
| LD clumping                   | `03_instrument_builder.R`      |
| Harmonisation                 | `04_outcome_harmonise.R`       |
| Sample overlap assessment     | `05_sample_overlap_check.R`    |
| MR estimation                 | `07_mr_estimation.R`           |
| pQTL validation               | `08_pqtl_validation.R`         |
| Colocalisation analysis       | `09_colocalisation.R`          |
| Benchmark classification      | `10_classification.R`          |
| Reporting tables              | `11_reporting_tables.R`        |
| Figures                       | `12_reporting_figures.R`       |
| Triangulation                 | `13_triangulation.R`           |
| Additional visualisation      | `14_additional_figures.R`      |
| Failure mode diagnostics      | `15_failure_mode_diagnostic.R` |
| Post-transcriptional analysis | `16_failure_mode_cat6.R`       |

---

## Output Structure

Each run generates a self-contained timestamped directory:

```
runs/
  <YYYYMMDD_HHMMSS>/
    results/
      *.csv
      figures/
    artifacts/
      config_snapshot.yaml
      session_info.txt
    pipeline_log.txt
    pipeline_errors.txt
```

---

## Failure Mode Diagnostic Framework

The pipeline includes a structured post-MR diagnostic system that systematically attributes null or discordant MR findings to specific methodological sources. Seven non-mutually exclusive diagnostic categories are evaluated for all *Not supported* and *Non-estimable* targets:

| Category | Domain |
|---|---|
| Category 1 | MHC/LD complexity |
| Category 2 | Structural variant / genomic cluster / paralog |
| Category 3 | pQTL panel absence |
| Category 4 | Tissue/cell-type specificity |
| Category 5 | Pharmacological mechanism mismatch |
| Category 6 | Post-transcriptional discordance (sQTL / eQTL–pQTL) |
| Category 7 | Instrument absence / statistical non-recovery |

Diagnostic flags are additive and confidence-tiered (High / Moderate / Low). They annotate each target without altering the primary benchmark classification.

Full category definitions, data sources, and decision rules are specified in `configs/failure_mode_diagnostic.yaml` and implemented in `R/15_failure_mode_diagnostic.R` and `R/16_failure_mode_cat6.R`.

---

## Dependencies

```r
install.packages(c(
  "dplyr", "tidyr", "yaml", "fs", "glue",
  "ggplot2", "data.table", "coloc", "susieR"
))

# TwoSampleMR from GitHub:
remotes::install_github("MRCIEU/TwoSampleMR")
```

| Package | Version (used) | Role |
|---|---|---|
| TwoSampleMR | 0.7.0 | Instrument formatting, harmonisation, MR estimation |
| coloc | 5.2.3 | Colocalisation analysis |
| susieR | 0.14.2 | SuSiE fine-mapping (triggered by protocol rules) |
| dplyr / tidyr | — | Data wrangling |
| yaml | — | Config file parsing |
| fs | — | File system utilities |
| glue | — | String interpolation in logging |
| ggplot2 | — | Figures |
| data.table | — | Sample overlap checks, large-file reading |
| MRInstruments | — | ARIES mQTL source data (**optional** — pre-extracted CSV included in `cache/`) |

---

## Reproducibility Guarantees

- Fully pre-specified protocol (YAML-locked before analysis)
- Config-driven execution (no hardcoded parameters)
- No manual intervention during analysis
- Automatic session capture and config snapshot per run
- Deterministic pipeline structure with attrition tracking

---

## Author

Mihye Kwon

---

## Citation

If you use this pipeline, please cite:

> Kwon M et al. (2026). btsDMARDs MR Pipeline (Version 1). Zenodo. https://doi.org/10.5281/zenodo.20051419

Manuscript citation to be added upon publication.

---

## License

MIT License — see [LICENSE](LICENSE)
