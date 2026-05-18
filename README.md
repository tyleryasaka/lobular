# lobular
Inference of hepatocyte zonation in a single cell (+/- spatial) dataset

The hepatic lobule demonstrates a gradient from the portal triad to the central vein, with hepatocytes (and other cell types) engaging in distinct gene expression and function along this gradient. These are broadly classified into zones 1, 2, and 3 based on their proximity to the portal vein (zone 1) or central vein (zone 3), with zone 2 denoting the intermediate region. (I have authored a [review of liver zonation](https://pubmed.ncbi.nlm.nih.gov/41082401/).)

The `lobular` R package evenly divides cells in a baseline sample into zones 1, 2, and 3 based on their expression of zonated genes. This is based on the spatial arrangement of genes along an annotated porto-central axis. This uses my own annotations which I performed on public data from the [Liver Cell Atlas](https://www.livercellatlas.org/). The majority of zonation captured is presumed to be hepatocyte derived, but this data is not hepatocyte-specific. There are numerous studies demonstrating that other cell types are zonated as well (e.g. endothelial cells, hepatic stellate cells, and macrophages; more detail on this in my review).

If you have any questions or comments, or any issues using the package, please feel free to reach out: `typublic@pm.me`

## Installation

```r
install.packages("devtools")
devtools::install_github("tyleryasaka/lobular")
```

## Usage

### Obtaining zonation in Seurat v5

```r
library(lobular)

# Train the model on a baseline sample
seurat_baseline <- subset(seurat_obj, subset = sample == 'my_baseline_sample_id')
zonation_obj <- trainModel(Seurat::GetAssayData(seurat_baseline, layer = 'data'), species = 'human')

# Apply the trained model to your full dataset (or any other sample)
zonation_obj <- applyModel(Seurat::GetAssayData(seurat_obj, layer = 'data'), zonation_obj)

# Obtain discrete zonation bins, normalized to the baseline
zonation_assignments <- getZone(zonation_obj)
seurat_obj <- AddMetaData(seurat_obj, zonation_assignments, col.name = 'zone')

# Obtain a continuous zonation gradient (values between 1 and 3), normalized to the baseline
zonation_gradient <- getZonationGradient(zonation_obj)
seurat_obj <- AddMetaData(seurat_obj, zonation_gradient, col.name = 'zonation')
```

## Functions

### Model setup

#### `trainModel()`

Train the model on a baseline liver sample.

**Usage:**
```r
trainModel(mtx, coords = NULL, species = 'human', regularization = 1, filter = 0, verbose = FALSE)
```

**Arguments:**
- `mtx`: Gene expression matrix (*log-normalized*) with genes as rows.
- `coords`: (Optional) For spatial data, calibrates to the dimensions of the baseline sample. Coordinate matrix with samples as rows, and columns `x` and `y`.
- `species`: (Optional) Species to use; defaults to `'human'`. Supports `'mouse'` and `'human'`.
- `regularization`: (Optional) Numeric regularization strength used during model fitting (default 1).
- `filter`: (Optional) Minimum absolute correlation threshold for genes to be retained in the model (default 0, i.e. no filtering).
- `verbose`: (Optional) If `TRUE`, prints training diagnostics (default `FALSE`).

**Returns:**
- A `ZonationObject` with calibrated baseline zonation.

---

#### `applyModel()`

Attach a (new) gene expression matrix to a trained `ZonationObject`. Almost every other function in the package operates on the resulting object.

**Usage:**
```r
applyModel(mtx, zone_obj)
```

**Arguments:**
- `mtx`: Gene expression matrix (*log-normalized*) with genes as rows.
- `zone_obj`: Trained `ZonationObject` (output of `trainModel()`).

**Returns:**
- A `ZonationObject` with `mtx` attached, ready for downstream queries and plotting.

---

### Extracting zonation

#### `getZone()`

Return the discrete zonation bin for each cell/spot.

**Usage:**
```r
getZone(zone_obj)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied (see `applyModel()`).

**Returns:**
- Factor of zonation assignments (`Zone_1`, `Zone_2`, `Zone_3`), one per cell/spot.

---

#### `getZonationGradient()`

Return the continuous zonation value (between 1 and 3) for each cell/spot.

**Usage:**
```r
getZonationGradient(zone_obj)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.

**Returns:**
- Numeric vector of zonation values on `[1, 3]`, one per cell/spot.

---

#### `getZonation2d()`

Return the 2-dimensional zonation representation (independent Zone 1 and Zone 3 scores) for each cell/spot.

**Usage:**
```r
getZonation2d(zone_obj)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.

**Returns:**
- A data frame with columns:
  - `ZONE_1`: Zone 1 (portal) score, scaled to `[0, 1]`.
  - `ZONE_3`: Zone 3 (central) score, scaled to `[0, 1]`.
  - `zonation`: The 1D zonation gradient value.
  - `zone`: Discrete zone assignment.

---

#### `getGeneZonation()`

Return per-gene zonation factors for each zone.

**Usage:**
```r
getGeneZonation(zone_obj)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject`.

**Returns:**
- A data frame with columns `zone_1`, `zone_2`, `zone_3` giving the per-gene zonation factor for each zone.

---

### Plotting zonation

#### `plotZonation2d()`

Scatter plot of cells/spots in the 2D Zone 1 × Zone 3 space, colored by discrete zone.

**Usage:**
```r
plotZonation2d(zone_obj, point_size = 1)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `point_size`: (Optional) Point size (default 1).

**Returns:**
- A `ggplot` object.

---

#### `plotZonation2d_2()`

Same as `plotZonation2d()`, but colored by Zone 2 score. **Mouse only** — there are no human Zone 2 reference genes.

**Usage:**
```r
plotZonation2d_2(zone_obj, point_size = 1)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied (`species = 'mouse'`).
- `point_size`: (Optional) Point size (default 1).

**Returns:**
- A `ggplot` object.

---

#### `plotZonation2dGene()`

Same as `plotZonation2d()`, but colored by the expression of a specified gene.

**Usage:**
```r
plotZonation2dGene(zone_obj, gene, point_size = 1)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `gene`: Name of the gene to plot (must be a row in the applied matrix).
- `point_size`: (Optional) Point size (default 1).

**Returns:**
- A `ggplot` object.

---

#### `plotRegression()`

Scatter and cubic polynomial fit of a gene's expression along the zonation axis.

**Usage:**
```r
plotRegression(zone_obj, gene)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `gene`: Name of the gene to plot (must be a row in the applied matrix).

**Returns:**
- A `ggplot` object.

---

#### `plotPolarity()`

2D density (heatmap) over Zone 1 × Zone 3 with the inferred polarity score annotated. Useful for visualizing how strongly anti-correlated the zone 1 and zone 3 axes are in a given sample.

**Usage:**
```r
plotPolarity(zone_obj)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.

**Returns:**
- A `ggplot` object.

---

#### `plotVirtualLobule()`

Plot an idealized hexagonal lobule, with each pixel placed by its rank along the corner-to-center axis (0 = nearest corner / Zone 1, 1 = center / Zone 3) and colored by the matching quantile of the inferred zonation gradient.

**Usage:**
```r
plotVirtualLobule(zone_obj,
                  resolution = 100,
                  palette = "ggthemes::Classic Red-Blue",
                  reverse_palette = TRUE,
                  pointy_top = FALSE,
                  show_legend = TRUE,
                  seed = NULL)
```

**Arguments:**
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `resolution`: (Optional) Integer pixel resolution along the x-axis (default 100).
- `palette`: (Optional) A `paletteer` continuous palette name (default `"ggthemes::Classic Red-Blue"`).
- `reverse_palette`: (Optional) If `TRUE`, reverse the palette direction (default `TRUE`).
- `pointy_top`: (Optional) If `TRUE`, the hexagon has a vertex at the top; otherwise flat-topped (default `FALSE`).
- `show_legend`: (Optional) If `TRUE`, show the color legend (default `TRUE`).
- `seed`: (Optional) Integer seed for reproducibility.

**Returns:**
- A `ggplot` object.

---

#### `plotZonationRidge()`

Overlay the inferred zonation distributions of one or more samples as outlined density curves.

**Usage:**
```r
plotZonationRidge(zone_objs, palette = 'grDevices::rainbow', line_width = 1, adjust = 3)
```

**Arguments:**
- `zone_objs`: A single `ZonationObject`, or a (preferably named) list of them. List names are used as legend labels.
- `palette`: (Optional) A `paletteer` palette name for the outline colors (default `'grDevices::rainbow'`).
- `line_width`: (Optional) Line width for density outlines (default 1).
- `adjust`: (Optional) Bandwidth adjustment passed to `geom_density` (default 3).

**Returns:**
- A `ggplot` object.

---

### Differential zonation

#### `plotZonationDiff()`

Plot the per-gene differential expression of zonated genes between two samples, with bar widths weighted by each gene's zonation factor.

**Usage:**
```r
plotZonationDiff(mtx_1, mtx_2, zone_obj, zone, threshold = 0.1, font_size = 9)
```

**Arguments:**
- `mtx_1`: Gene expression matrix (*log-normalized*) for sample 1 (reference). Genes as rows.
- `mtx_2`: Gene expression matrix (*log-normalized*) for sample 2 (comparison). Genes as rows.
- `zone_obj`: Calibrated `ZonationObject`.
- `zone`: Numeric zone to plot. Must be `1` or `3` for human, or `1`, `2`, or `3` for mouse.
- `threshold`: (Optional) Minimum zonation factor for a gene to be included (default 0.1).
- `font_size`: (Optional) X-axis label font size (default 9).

**Returns:**
- A `ggplot` object showing per-gene expression change (sample 2 − sample 1) for the requested zone.

---

### Spatial zonation

These functions infer zonation from gene expression and then propagate it across a tissue section via spatial interpolation. This is particularly useful when you want to assign zonation to non-parenchymal cells based on their proximity to zonated hepatocytes rather than their own expression.

#### `getZoneSpatial()`

Infer discrete zonation for all cells/spots based on spatial interpolation.

**Usage:**
```r
getZoneSpatial(coords, zone_obj, resolution = 1, use_for_inference = NULL)
```

**Arguments:**
- `coords`: Coordinate matrix with samples as rows, and columns `x` and `y`. Rownames should match the colnames of the matrix applied to `zone_obj`.
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `resolution`: (Optional) Interpolation granularity (default 1; higher = finer).
- `use_for_inference`: (Optional) A vector of sample names which should be used for zonation inference (recommended: hepatocytes only, if annotation is available). If not provided, all samples are used.

**Returns:**
- Factor of discrete zonation assignments for all samples.

---

#### `plotZoneSpatial()`

Plot the 2D interpolated zones as a filled contour map.

**Usage:**
```r
plotZoneSpatial(coords, zone_obj, resolution = 1, use_for_inference = NULL)
```

**Arguments:**
- `coords`: Coordinate matrix with columns `x` and `y`.
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `resolution`: (Optional) Interpolation granularity (default 1).
- `use_for_inference`: (Optional) Subset of samples for inference.

**Returns:**
- A `ggplot` object with filled contour zones.

---

#### `plotZoneSpatialContours()`

Plot 2D zonation with contour outlines overlaid on the original points.

**Usage:**
```r
plotZoneSpatialContours(coords, zone_obj, resolution = 1, point_size = 1, line_width = 2, plot_options = NULL, use_for_inference = NULL)
```

**Arguments:**
- `coords`: Coordinate matrix with columns `x` and `y`.
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `resolution`: (Optional) Interpolation granularity (default 1).
- `point_size`: (Optional) Point size (default 1).
- `line_width`: (Optional) Contour line width (default 2).
- `plot_options`: (Optional) Reserved for additional plot customization.
- `use_for_inference`: (Optional) Subset of samples for inference.

**Returns:**
- A `ggplot` object with points colored by zonation and contour outlines.

---

#### `plotZoneSpatialCustom()`

Plot a custom variable spatially, with zonation contour outlines overlaid.

**Usage:**
```r
plotZoneSpatialCustom(meta, col_name, zone_obj, resolution = 1, point_size = 1, use_for_inference = NULL)
```

**Arguments:**
- `meta`: Metadata data frame with samples as rows, and columns `x`, `y`, and the column named by `col_name`. Rownames should match the colnames of the matrix applied to `zone_obj`.
- `col_name`: Name of the column in `meta` to color points by.
- `zone_obj`: Calibrated `ZonationObject` with a matrix applied.
- `resolution`: (Optional) Interpolation granularity (default 1).
- `point_size`: (Optional) Point size (default 1).
- `use_for_inference`: (Optional) Subset of samples for inference.

**Returns:**
- A `ggplot` object with points colored by `col_name` and contour outlines.
