# In situ scanner-based imaging to monitor edaphic activity

This repository supports the development of a **soil biological activity monitoring method** based on buried flatbed scanners. CIS-type flatbed scanners are installed vertically in the soil and continuously capture high-resolution images of the soil matrix *in situ*, enabling non-destructive, long-term observation of soil organisms — including invertebrates and root structures — at fine spatial and temporal resolutions.

Coupled with automated image capture (via Raspberry Pi) and a suite of image processing and analysis tools, this method opens new possibilities for studying edaphic biodiversity across diverse spatio-temporal scales ([Belaud et al. 2024](https://doi.org/10.1007/s00374-024-01851-8)).

![](image/scan-device.png)

------------------------------------------------------------------------

# Overview

The methodology spans two main components:

| Component | Description |
|---------------------------------|---------------------------------------|
| **Image Acquisition** | Scanner setup, automated capture, resolution/frequency calibration |
| **Image Analysis** | Invertebrate & root detection/classification via computer vision and deep learning |

------------------------------------------------------------------------

# Repository Structure

```         
soil_scan/
├── image/                        # Repository figures and device photos
├── images_acquisition/           # Acquisition calibration code and data
│   ├── 00_run_all.R              # Full pipeline runner (start here)
│   ├── 01_functions.R            # Helper functions (parsing, run IDs)
│   ├── 02_prepare_fit_models.R   # Data preparation and model fitting (GLMMs)
│   ├── 03_tables_figures.R       # Output generation (tables and figures)
│   ├── container-sensors.R       # Sensor container utilities
│   ├── data/                     # Raw calibration data
│   │   ├── database_calibration.csv   # Labelled invertebrate detections
│   │   └── BACs_sensors.dat           # Environmental sensor readings
│   ├── Tables/                   # Generated summary tables (CSV)
│   └── Figures/                  # Generated figures (PNG, PDF)
└── images_analysis/              # Image analysis pipeline related study
    ├── Use-case.qmd              # Reproducible ecological use-case analysis
    ├── pipeline_evaluation.qmd   # Pipeline performance evaluation
    ├── DESCENT_data/             # DESCENT project monitoring data
    ├── DIAMS_data/               # DIAMS project monitoring data
    ├── TRUF_data/                # TRUF project monitoring data
    ├── css/                      # Stylesheet for HTML reports
    └── output/                   # Generated figures
```

------------------------------------------------------------------------

## Image acquisition

> **Paper:** *"In situ scanner-based imaging to monitor soil mesofauna activity: trade-offs between spatial resolution and temporal frequency"*

📁 Code: [`images_acquisition/`](images_acquisition/)

A core challenge in deploying buried scanner systems is that image-derived activity metrics are shaped by *how images are acquired*, not just by biological activity itself. This study experimentally calibrates the effects of two key acquisition parameters on the detectability of soil mesofauna under **controlled greenhouse conditions**:

-   **Spatial resolution** — acts as a size-selective filter, determining which organisms and size classes enter the observation window
-   **Capture frequency** — modulates detection magnitude through temporal integration, without redefining the detectable size spectrum

Three taxonomic groups are studied: **Acari**, **Collembola**, and **Enchytraeidae**.

These findings provide a mechanistic framework for how protocol choices translate into observable signals, and offer practical guidance for standardising acquisition settings — a prerequisite for robust ecological inference and cross-study comparability.

### Running the Acquisition Analysis

``` r
# From the project root, run the full pipeline:
source("images_acquisition/00_run_all.R")
```

Scripts must be run in order (handled automatically by `00_run_all.R`):

1.  `01_functions.R` — Loads helper functions
2.  `02_prepare_fit_models.R` — Prepares data and fits GLMMs
3.  `03_tables_figures.R` — Generates all tables and figures

------------------------------------------------------------------------

## Image Analysis

> **Paper:** *"Automated processing for in situ soil monitoring: from raw images to population dynamics"*

📁 Code: [`images_analysis/`](images_analysis/) · 📄 [View the analysis online](https://htmlpreview.github.io/?https://github.com/emmabelaud/soil_scan/blob/24b6554309cb9a2cf67c79a98833a7e2e406e59e/images_analysis/Use-case.html)

This component develops an **end-to-end processing pipeline** combining traditional computer vision with state-of-the-art deep learning, designed to operate under the constrained annotation budgets typical of ecological studies.

Key contributions:

-   **Image differencing** reframes the detection problem (complicated by very low signal-to-noise ratios in soil imagery) into a simpler classification task
-   **Fine-tuned foundation models** trained on limited labelled data deliver reliable population counts across taxa
-   A **labelled dataset** of \~600 soil images and \>8,000 annotated invertebrates, spanning Acari, Collembola, Enchytraeidae, and others
-   Application to **4 scanners over 3 months** of continuous monitoring, demonstrating longitudinal biodiversity analyses at previously unreachable scales

### Running the Analysis Reports

Requires [Quarto](https://quarto.org/) to be installed.

``` bash
# Render the ecological use-case report
quarto render images_analysis/Use-case.qmd

# Render the pipeline evaluation report
quarto render images_analysis/pipeline_evaluation.qmd
```

------------------------------------------------------------------------

## Getting Started

1.  **Clone the repository**

``` bash
git clone https://github.com/emmabelaud/soil_scan.git
cd soil_scan
```

2.  **Open the R project**

Open `soil_scan.Rproj` in RStudio to set the working directory correctly.

3.  **Run the acquisition analysis**

``` r
source("images_acquisition/00_run_all.R")
```

4.  **Render the analysis reports**

``` bash
quarto render images_analysis/Use-case.qmd
```

------------------------------------------------------------------------

## Data

| File | Location | Description |
|-----------------|------------------------|-------------------------------|
| `database_calibration.csv` | `images_acquisition/data/` | Labelled invertebrate detections with resolution/frequency metadata |
| `BACs_sensors.dat` | `images_acquisition/data/` | Environmental sensor readings from scanner containers |
| `SLM_2023_combined_pipeline.csv` | `images_analysis/DESCENT_data/` | DESCENT project pipeline output |
| `combined_data.csv` | `images_analysis/DIAMS_data/` | DIAMS project pipeline output |
| `TRUF_data_combined_pipeline.csv` | `images_analysis/TRUF_data/` | TRUF project pipeline output |
