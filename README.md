# lobular
Inference of hepatocyte zonation in a single cell (+/- spatial) dataset

The hepatic lobule demonstrates a gradient from the portal triad to the central vein, with hepatocytes (and other cell types) engaging in distinct gene expression and function along this gradient. These are broadly classified into zones 1, 2, and 3 based on their proximity to the portal vein (zone 1) or central vein (zone 3), with zone 2 denoting the intermediate region.

The `lobular` R package evenly divides cells in a baseline sample into zones 1, 2, and 3 based on their expression of zonated genes, as previously published ([Xu et al, *Nature Genetics* 2024](https://pubmed.ncbi.nlm.nih.gov/38627598/)). The zonation assignments are based on spatial transcriptomics data across the entire lobule. The majority of zonation captured is presumed to be hepatocyte derived, but this data is not hepatocyte-specific. There are numerous studies demonstrating that other cell types are zonated as well (e.g. endothelial cells, hepatic stellate cells, and macrophages).

## Installation

```r
install.packages("devtools")
devtools::install_github("tyleryasaka/lobular")
```

## Loading

```r
library(lobular)
```

## Usage

1. Call `setBaseline(mtx, species)` with your gene expression matrix and desired species to calibrate the model and obtain a `ZonationObject`.

2. Use the calibrated `ZonationObject` and your gene expression matrix as input to `getZone(mtx, zone_obj)` to obtain zonation assignments for each sample in your matrix.

## Examples

### Obtaining zonation in Seurat v5

```r
seurat_baseline <- subset(seurat_obj, subset = sample == 'my_baseline_sample_id')
zonation_obj <- setBaseline(Seurat::GetAssayData(seurat_baseline, slot = 'data'))
zonation_assignments <- getZone(Seurat::GetAssayData(seurat_obj, slot = 'data'), zonation_obj)
seurat_obj <- AddMetaData(seurat_obj, zonation_assignments, col.name = 'zone')
```

## Functions

### `setBaseline(mtx, species = 'human')`

Calibrate the model to baseline liver zonation.

#### Parameters

- `mtx`: Gene expression matrix with genes as rows
- `species`: Species to use, defaults to human. Currently supports 'mouse' and 'human'.

#### Return Value

A `ZonationObject` with calibrated baseline zonation.

### `getZone(mtx, zone_obj)`

Apply the model to new values, returning the zonation.

#### Parameters

- `mtx`: Gene expression matrix with genes as rows
- `zone_obj`: Calibrated Zonation Object

#### Return Value

A vector of zonation assignments. (1, 2, and 3)
