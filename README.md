# btsDMARDs MR Pipeline (Version 1)

## Overview

This repository provides a fully reproducible Mendelian randomisation (MR) pipeline for evaluating drug target validity using multi-layer genetic evidence, including eQTL, pQTL, colocalisation, and triangulation frameworks.

The pipeline is **protocol-driven, pre-specified, and designed for publication-grade reproducibility**, enabling transparent and auditable causal inference in genetic epidemiology.

---

## Scientific Objective

The primary objective of this pipeline is to:

* Evaluate causal effects of genetically proxied drug targets
* Integrate multi-omics evidence (eQTL, pQTL, colocalisation)
* Identify mechanistic failure modes in MR inference
* Provide a structured triangulation framework for drug target validation

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

---

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

  R/
    00_utils.R
    01_protocol_checks.R
    02_input_readers.R
    03_instrument_builder.R
    04_outcome_harmonise.R
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
    protocol_v1.yaml
    run_full_study_v1.yaml
    failure_mode_diagnostic.yaml
    targets.yaml
    pqtl_sources.yaml
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

## How to Run

```r
source("main_fullstudy.R")

cfg <- read_protocol_bundle(
  "configs/protocol_v1.yaml",
  "configs/run_full_study_v1.yaml"
)
```

---

## External dependencies

This pipeline requires:

- PLINK 1.9
- 1000 Genomes reference panels (EUR/EAS)

Set environment variables:

- PLINK_BIN
- LD_REF_EUR
- LD_REF_EAS
- INPUT_DATA_PATH

---

## Data Input

This repository does not include raw data.

To run the pipeline:

### Option 1: Use environment variable

```r
Sys.setenv(INPUT_DATA_PATH = "/path/to/data.csv")
```

### Option 2: Place file locally

```
data/raw/input.csv
```

## Exposure Data Processing

Exposure data were derived from the eQTLGen cis-eQTL summary statistics.

Effect sizes (beta) and standard errors (SE) were computed from Z-scores using:

- beta = Z / sqrt(N)
- SE = 1 / sqrt(N)

where Z is the reported Z-score and N is the sample size.

Allele frequencies were obtained from the eQTLGen SNP allele frequency file:

`2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz`

Available at:  
https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz

The effect allele frequency (EAF) corresponds to `AlleleB_all` and was aligned to the assessed allele used in the eQTL summary statistics.

The processed exposure dataset used in this pipeline:

`exposure_eqtlgen_targets18_betaSE.txt`

All transformations were performed prior to analysis; no additional transformations are applied within the pipeline.


---

## Data Availability

External datasets are not distributed with this repository due to licensing restrictions. Users must obtain the datasets independently and configure paths via environment variables.

---

## Output

Output files will be written to the `runs/` directory, which is created automatically during execution.

Each run generates a timestamped folder containing results, logs, and reproducibility artifacts.

---

## Output Structure

```
runs/
  eqtl_abstract_<timestamp>/
    results/
      *.csv
      figures/
    artifacts/
      config_snapshot.yaml
      session_info.txt
    pipeline_log.txt
```

---

## Failure Mode Diagnostic Framework

The pipeline includes a structured post-MR diagnostic system:

* Category 7: Instrument absence (IV failure cascade)
* Category 6: Post-transcriptional discordance
* Category 5: Mechanism mismatch
* Category 3: pQTL panel absence
* Category 1–2: LD structure / genomic complexity

This framework enables **systematic interpretation of null or discordant MR findings**.

---

## Dependencies

* R (≥ 4.0)
* dplyr
* TwoSampleMR
* coloc
* yaml
* fs
* ggplot2

---

## Reproducibility Guarantees

* Fully pre-specified protocol
* Config-driven execution
* No manual intervention during analysis
* Automatic logging and session capture
* Deterministic pipeline structure

---

## Author

Mihye Kwon

---

## Citation

(To be added upon publication)

---

## License

(To be specified)
