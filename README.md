<div align="center">

# Online Ant Trade

### Risk Assessment of Ant Imports to Norway

*Species distribution modeling pipeline supporting VKM's risk assessment of importing ant species for private hobby keeping in Norway.*

![Status](https://img.shields.io/badge/status-active-brightgreen)
![Language](https://img.shields.io/badge/language-R-276DC3?logo=r)
![Method](https://img.shields.io/badge/method-MaxEnt%20%2F%20SDMtune-orange)
![Climate](https://img.shields.io/badge/climate-WorldClim%20SSP585-blue)
![Commissioned by](https://img.shields.io/badge/commissioned%20by-Milj%C3%B8direktoratet-2E7D32)

</div>

---

## Table of Contents

- [Background](#background)
- [What this repository does](#what-this-repository-does)
- [Repository layout](#repository-layout)
- [Pipeline highlights](#pipeline-highlights)
- [Getting started](#getting-started)
- [Commissioning context](#commissioning-context)
- [References](#references)
- [License](#license)

---

## Background

> The Norwegian Environment Agency (**Miljødirektoratet**) commissioned VKM (ref. **2025/15126**, 16 October 2025) to assess the risk of adverse effects on biological diversity from the import of various ant species for **private hobby keeping** or trade to hobby keepers.

Miljødirektoratet has received applications for permits to import close to **100 ant species**, spanning tropical to boreal distributions, including genera such as *Camponotus*, *Crematogaster*, *Messor*, and *Myrmecia*. The international hobby market offers an even larger pool of species than has been formally requested.

Under the Norwegian Regulation on Alien Organisms (*forskrift om fremmede organismer* § 6), Miljødirektoratet must evaluate whether each requested species poses a risk of adverse effects on biological diversity, including whether it can **survive and establish in Norway**. There is a general lack of knowledge about the establishment potential of these ant species under Norwegian conditions and the possible consequences for native biodiversity.

VKM's assessment of establishment ability and possible risk to biological diversity forms part of the knowledge base that Miljødirektoratet uses to decide whether species can be permitted, banned, or exempted from permit requirements under the regulation.

---

## What this repository does

A reproducible R / Python pipeline that:

| Stage | Description |
|---|---|
| **Data acquisition** | Occurrences from GBIF, BISON, ALA, iNaturalist, iDigBio, and the GABI AntMaps database |
| **Cleaning** | `CoordinateCleaner` filters + spatial thinning to one record per raster cell |
| **Predictors** | WorldClim BIO, soil-bioclim (SBIO), and additional environmental layers |
| **Variable filtering** | Keeps all BIO; removes SBIO/ENV correlated with BIO (Spearman \|r\| > 0.7) |
| **Modeling** | MaxEnt via `SDMtune` with cross-validated hyperparameter tuning |
| **Projection** | Current and future climate (SSP585, 2021–2040), clamped + unclamped |
| **Diagnostics** | Cross-validation SD, MESS extrapolation, clamping-effect maps |
| **Outputs** | Maps, response curves, ROC curves, variable importance, per-species summaries |

The output supports VKM's evaluation of which ant species can plausibly establish in Norway under current and projected climates.

---

## Repository layout

```
.
├── 4_scripts/                  R / Python pipeline scripts
│   └── Gamle skript/           Older script versions kept for reference
├── presentations/              Project presentations and templates
├── CHANGELOG.md                Pipeline feature / fix history
├── README.md                   This file
└── .gitignore                  Excludes data, outputs, PDFs, Claude files, etc.
```

> **Note** — Data folders (`1_raw_data/`, `2_processed_data/`, `5_outputs/`, `species/`, `pdf/`, `reports/`, …) and binary outputs (`*.tif`, `*.pdf`) are intentionally **not** version-controlled. See `.gitignore`.

---

## Pipeline highlights

The production SDM script is **`4_scripts/10_updated_sdmtune_loop_merged.R`**. Full feature history lives in [`CHANGELOG.md`](CHANGELOG.md).

- **Clamping** — clamped *and* unclamped predictions for current and future, with diagnostic difference maps
- **Smart variable filtering** — keep all BIO climate variables; drop SBIO/ENV that are redundant with BIO; let `varSel` choose among the remainder
- **GABI AntMaps integration** — local GABI occurrences merged with online sources, matched by species name *and* GBIF acceptedKey
- **Ensemble modeling** — standalone example combining Maxnet, Random Forest, and BRT via AUC-weighted mean
- **Robust background sampling** — biogeographically-informed circular buffers with automatic fallback to random global background
- **Timing instrumentation** — per-species and overall pipeline timing recorded in summary tables

---

## Getting started

### Prerequisites

R ≥ 4.2 with the following packages:

```r
install.packages(c(
  "spocc", "rgbif", "tidyverse", "CoordinateCleaner",
  "SDMtune", "terra", "sf", "ggplot2", "rio", "rnaturalearth",
  "raster", "dismo", "patchwork", "viridis", "tidyterra"
))
```

### Run the production pipeline

```r
source("./4_scripts/10_updated_sdmtune_loop_merged.R")
```

Configure species list, minimum records, CV folds, and buffer distance at the top of the script.

---

## Commissioning context

| | |
|---|---|
| **Commissioning agency** | Miljødirektoratet — Norwegian Environment Agency |
| **Performing institution** | VKM — Vitenskapskomiteen for mat og miljø |
| **Reference** | 2025/15126 |
| **Date of commission** | 16 October 2025 |
| **Legal basis** | Forskrift om fremmede organismer § 6 |
| **Case officer** | Emma Erlingsson (Miljødirektoratet) |

---

## References

Key methodological references underpinning the pipeline:

- **Elith, J. et al. (2011)** — *A statistical explanation of MaxEnt for ecologists.* Diversity and Distributions 17(1): 43–57.
- **Merow, C., Smith, M.J., & Silander, J.A. (2013)** — *A practical guide to MaxEnt for modeling species' distributions.* Ecography 36(10): 1058–1069.
- **Radosavljevic, A. & Anderson, R.P. (2014)** — *Making better Maxent models of species distributions: complexity, overfitting and evaluation.* Journal of Biogeography 41(4): 629–643.
- **Vignali, S. et al. (2020)** — *SDMtune: An R package to tune and evaluate species distribution models.* Ecology and Evolution 10(20): 11488–11506.
- **Zizka, A. et al. (2019)** — *CoordinateCleaner: Standardized cleaning of occurrence records from biological collection databases.* Methods in Ecology and Evolution 10(5): 744–751.

---

## License

To be determined.

<div align="center">

---

*Maintained by Daniel Flo (VKM) — for questions, please open an issue.*

</div>
