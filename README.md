# Online Ant Trade — Risk Assessment of Ant Imports to Norway

Species distribution modeling and data pipeline supporting VKM's (Norwegian Scientific Committee for Food and Environment) risk assessment of importing ant species for private hobby keeping in Norway.

## Background

The Norwegian Environment Agency (Miljødirektoratet) commissioned VKM (ref. 2025/15126, 16 October 2025) to assess the risk of adverse effects on biological diversity from the import of various ant species for private hobby keeping or trade to hobby keepers.

Miljødirektoratet has received applications for permits to import a large number of ant species — close to 100 species so far — spanning tropical to boreal distributions, including genera such as *Camponotus*, *Crematogaster*, *Messor*, and *Myrmecia*. The international hobby market offers an even larger pool of species than has been formally requested.

Under the Norwegian Regulation on Alien Organisms (*forskrift om fremmede organismer* § 6), Miljødirektoratet must evaluate whether each requested species poses a risk of adverse effects on biological diversity, including whether it can survive and establish in Norway. There is a general lack of knowledge about the establishment potential of these ant species under Norwegian conditions and the possible consequences for native biodiversity.

VKM's assessment of establishment ability and possible risk to biological diversity forms part of the knowledge base that Miljødirektoratet uses to decide whether species can be permitted, banned, or exempted from permit requirements under the regulation.

## What this repository contains

R (and a few Python) scripts implementing a reproducible pipeline that:

- Acquires occurrence data for ant species from GBIF, BISON, ALA, iNaturalist, iDigBio, and the GABI AntMaps database
- Cleans and thins occurrence records (CoordinateCleaner + spatial thinning)
- Prepares bioclimatic, soil-bioclim, and other environmental predictor layers
- Filters correlated environmental variables (keeps all BIO; removes redundant SBIO/ENV)
- Fits MaxEnt species distribution models with `SDMtune` (cross-validated, hyperparameter-tuned)
- Projects suitability to current and future climates (SSP585, 2021–2040), with both clamped and unclamped predictions
- Produces uncertainty (CV-fold SD), MESS extrapolation, and clamping-effect diagnostics
- Generates maps, response curves, ROC curves, variable importance plots, and per-species summary tables

The output is intended to support VKM's evaluation of which ant species can plausibly establish in Norway under current and projected climates.

## Repository layout

```
.
├── 4_scripts/                  R/Python pipeline scripts
│   └── Gamle skript/           Older script versions kept for reference
├── presentations/              Project presentations and templates
├── CHANGELOG.md                Pipeline feature/fix history
├── README.md                   This file
└── .gitignore                  Excludes data, outputs, PDFs, Claude files, etc.
```

Data folders (`1_raw_data/`, `2_processed_data/`, `3_metadata/`, `5_outputs/`, `species/`, `pdf/`, `reports/`, etc.) and binary outputs (`*.tif`, `*.pdf`) are intentionally **not** version-controlled — see `.gitignore`.

## Pipeline overview

The production SDM script is `4_scripts/10_updated_sdmtune_loop_merged.R`. See `CHANGELOG.md` for the consolidated history of features:

- **Clamping**: clamped + unclamped predictions for current and future, with diagnostic difference maps
- **Variable filtering**: keep all BIO variables, remove SBIO/ENV correlated with BIO (Spearman |r| > 0.7), then `varSel` selects among BIO
- **GABI AntMaps integration**: local GABI occurrences merged with online sources, matched by name and GBIF acceptedKey
- **Ensemble modeling**: standalone example combining Maxnet, Random Forest, and BRT via AUC-weighted mean
- **Timing instrumentation**: per-species and overall pipeline timing recorded in summary tables

## Required R packages

```r
install.packages(c(
  "spocc", "rgbif", "tidyverse", "CoordinateCleaner",
  "SDMtune", "terra", "sf", "ggplot2", "rio", "rnaturalearth",
  "raster", "dismo", "patchwork", "viridis", "tidyterra"
))
```

## Commissioning context

| | |
|---|---|
| **Commissioning agency** | Miljødirektoratet (Norwegian Environment Agency) |
| **Performing institution** | VKM — Vitenskapskomiteen for mat og miljø |
| **Reference** | 2025/15126 |
| **Date of commission** | 16 October 2025 |
| **Legal basis** | Forskrift om fremmede organismer § 6 |
| **Case officer** | Emma Erlingsson (Miljødirektoratet) |

## License

To be determined.
