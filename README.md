# lobular
Inference of hepatocyte zonation in a single cell (+/- spatial) dataset

The hepatic lobule demonstrates a gradient from the portal triad to the central vein, with hepatocytes (and other cell types) engaging in distinct gene expression and function along this gradient. These are broadly classified into zones 1, 2, and 3 based on their proximity to the portal vein (zone 1) or central vein (zone 3), with zone 2 denoting the intermediate region. (I have authored a [review of liver zonation](https://pubmed.ncbi.nlm.nih.gov/41082401/).)

The `lobular` R package evenly divides cells in a baseline sample into zones 1, 2, and 3 based on their expression of zonated genes. This is based on the spatial arrangement of genes along an annotated porto-central axis. For *mouse*, this uses previously published annotations ([Xu et al, *Nature Genetics* 2024](https://pubmed.ncbi.nlm.nih.gov/38627598/)). For *human*, this uses my own annotations which I performed on public data from the [Liver Cell Atlas](https://www.livercellatlas.org/). The majority of zonation captured is presumed to be hepatocyte derived, but this data is not hepatocyte-specific. There are numerous studies demonstrating that other cell types are zonated as well (e.g. endothelial cells, hepatic stellate cells, and macrophages; more detail on this in my review).

## Installation

```r
install.packages("devtools")
devtools::install_github("tyleryasaka/lobular")
```

## Usage

### Obtaining zonation in Seurat v5

```r
library(lobular)

# Normalize to a baseline sample
seurat_baseline <- subset(seurat_obj, subset = sample == 'my_baseline_sample_id')
zonation_obj <- setBaseline(Seurat::GetAssayData(seurat_baseline, layer = 'data'), species='human')

# Obtain discrete zonation bins for entire dataset, normalized to the baseline
zonation_assignments <- getZone(Seurat::GetAssayData(seurat_obj, layer = 'data'), zonation_obj)
seurat_obj <- AddMetaData(seurat_obj, zonation_assignments, col.name = 'zone')

# Obtain a continuous zonation gradient for entire dataset, normalized to the baseline
zonation_gradient <- getZonationGradient(Seurat::GetAssayData(seurat_obj, layer = 'data'), zonation_obj)
seurat_obj <- AddMetaData(seurat_obj, zonation_gradient, col.name = 'zonation')
```

## Functions

### `setBaseline()`

Calibrate the model to baseline liver zonation.

**Usage:**
```r
setBaseline(mtx, species = 'human')
```

**Arguments:**
- `mtx`: Gene expression matrix with genes as rows.
- `coords`: (Optional) For spatial data, calibrates to the dimensions of the baseline sample. Coordinate matrix with samples as rows, and columns `x` and `y`.
- `species`: (Optional) Species to use; defaults to `'human'`. Supports `'mouse'` and `'human'`.

**Returns:**
- A `ZonationObject` with calibrated baseline zonation.

---

### `getZone()`

Apply the model to new samples, returning discrete zonation bins.

**Usage:**
```r
getZone(mtx, zone_obj)
```

**Arguments:**
- `mtx`: Gene expression matrix with genes as rows.
- `zone_obj`: Calibrated `ZonationObject`.

**Returns:**
- Vector of zonation assignments (`Zone_1`, `Zone_2`, `Zone_3`) as a factor.

---

### `getZonationGradient()`

Apply the model to new samples, returning continuous zonation (0–1).

**Usage:**
```r
getZonationGradient(mtx, zone_obj)
```

**Arguments:**
- `mtx`: Gene expression matrix with genes as rows.
- `zone_obj`: Calibrated `ZonationObject`.

**Returns:**
- Vector of numeric zonation values (continuous) per sample.

---

### `getZoneSpatial()`

Infer zonation for all cells based on spatial interpolation.

**Usage:**
```r
getZoneSpatial(mtx, coords, zone_obj, use_for_inference = NULL)
```

**Arguments:**
- `mtx`: Gene expression matrix with genes as rows.
- `coords`: Coordinate matrix with columns `x` and `y`; rownames match `mtx` colnames.
- `zone_obj`: Calibrated `ZonationObject`.
- `resolution`: Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
- `use_for_inference`: Optional vector of sample names to use for zonation inference (e.g., hepatocytes).

**Returns:**
- Vector of discrete zonation assignments for all samples.

---

### `plotZoneSpatial()`

Plot interpolated zonation zones in 2D.

**Usage:**
```r
plotZoneSpatial(mtx, coords, zone_obj, use_for_inference = NULL)
```

**Arguments:**
- `mtx`: Gene expression matrix with genes as rows.
- `coords`: Coordinate matrix with columns `x` and `y`.
- `zone_obj`: Calibrated `ZonationObject`.
- `resolution`: Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
- `use_for_inference`: Optional subset of samples for inference.

**Returns:**
- A `ggplot` object with filled contour zones.

---

### `plotZoneSpatialContours()`

Plot 2D zonation with contour outlines and points.

**Usage:**
```r
plotZoneSpatialContours(mtx, coords, zone_obj, use_for_inference = NULL)
```

**Arguments:**
- `mtx`: Gene expression matrix with genes as rows.
- `coords`: Coordinate matrix with columns `x` and `y`.
- `zone_obj`: Calibrated `ZonationObject`.
- `resolution`: Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
- `point_size`: Optional numeric value for the ggplot point size (default 1)
- `use_for_inference`: Optional subset of samples for inference.

**Returns:**
- A `ggplot` object with points colored by zonation and contour outlines.

---

### `plotZoneSpatialCustom()`

Plot a custom variable with zonation contour outlines.

**Usage:**
```r
plotZoneSpatialCustom(mtx, meta, zone_obj, use_for_inference = NULL)
```

**Arguments:**
- `mtx`: Gene expression matrix with genes as rows.
- `meta`: Metadata matrix with samples as rows, and columns `x`, `y`, and `mycolname`, where `mycolname` is passed as `colname`. Rownames of coords should match colnames of mtx.
- `colname`: Name of custom column in `meta`
- `zone_obj`: Calibrated `ZonationObject`.
- `resolution`: Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
- `point_size`: Optional numeric value for the ggplot point size (default 1)
- `use_for_inference`: Optional subset of samples for inference.

**Returns:**
- A `ggplot` object with points colored by `label` and contour outlines.
