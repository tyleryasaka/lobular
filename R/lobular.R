ZONES = c('Zone_1', 'Zone_2', 'Zone_3')

#' Create a Zonation Object
#'
#' @param baseline Numeric vector of baseline zonation scores
#' @param species Character string indicating species ("human" or "mouse")
#' @param factors Named numeric vector of zonation factors for genes
#' @return An object of class \code{ZonationObject} containing:
#'   \item{baseline}{Baseline zonation scores}
#'   \item{species}{Species used for analysis}
#'   \item{factors}{Gene-level zonation factors}
#' @export
ZonationObject = function(baseline, baseline_pv, baseline_cv, species, factors, factors_pv, factors_cv, scale_factor) {
  # Input validation
  if (!is.numeric(baseline)) {
    stop("baseline must be numeric")
  }
  if (!species %in% c("human", "mouse")) {
    stop("species must be 'human' or 'mouse'")
  }
  if (!is.numeric(factors) || is.null(names(factors))) {
    stop("factors must be a named numeric vector")
  }
  if (!is.numeric(scale_factor) || is.null(scale_factor)) {
    stop("scale_factor must be a number")
  }

  # Create the object
  obj = list(
    baseline = baseline,
    baseline_pv = baseline_pv,
    baseline_cv = baseline_cv,
    species = species,
    factors = factors,
    factors_pv = factors_pv,
    factors_cv = factors_cv,
    scale_factor = scale_factor
  )

  class(obj) = "ZonationObject"
  return(obj)
}

#' Print method for ZonationObject
#' @param x A ZonationObject
#' @param ... Additional arguments (ignored)
#' @export
print.ZonationObject = function(x, ...) {
  cat("Zonation Object\n")
  cat("Species:", x$species, "\n")
  cat("Number of genes:", length(x$factors), "\n")
  cat("Baseline samples:", length(x$baseline), "\n")
  cat("Sale factor:", x$scale_factor, "\n")
  invisible(x)
}

#' Summary method for ZonationObject
#' @param object A ZonationObject
#' @param ... Additional arguments (ignored)
#' @export
summary.ZonationObject = function(object, ...) {
  cat("Zonation Object Summary\n")
  cat("=======================\n")
  cat("Species:", object$species, "\n")
  cat("Number of genes:", length(object$factors), "\n")
  cat("Baseline samples:", length(object$baseline), "\n")
  cat("\nBaseline distribution:\n")
  print(summary(object$baseline))
  cat("\nFactor range:", range(object$factors), "\n")
  cat("Sale factor:", object$scale_factor, "\n")
  invisible(object)
}

#' Convert mouse genes to human
#'
#' @param mouse_genes A vector of mouse gene names
#' @return A vector of human gene names
#' @noRd
mouseToHuman = function(mouse_genes) {
  mouse_genes_valid = mouse_genes[!is.na(mouse_genes)]
  conversion = orthogene::convert_orthologs(gene_df=mouse_genes_valid, method = 'gprofiler', gene_output='dict', input_species='mouse', output_species='human', non121_strategy='keep_popular', agg_fun = 'sum')
  sapply(mouse_genes, function(gene) ifelse(gene %in% mouse_genes_valid, conversion[gene], NA))
}

#' Obtain a weighted average of genes from a gene expression matrix
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param gene_weights A named vector of weights for each gene by gene name
#' @return A vector of weighted gene expression
#' @noRd
getGeneAvg = function(mtx, gene_weights) {
  genes = intersect(rownames(mtx), names(gene_weights))
  if (length(genes) == 0) {
    stop('No overlapping genes between Seurat object and gene_weights')
  }
  sub_expr = mtx[genes, , drop = FALSE]
  w = gene_weights[genes]
  scores = Matrix::colMeans(sub_expr * w)
  return(scores)
}

#' Transform a new vector to the distribution of an existing vector
#'
#' @param new_values Vector of new values
#' @param original_values Vector of original values
#' @return The transformed new vector
#' @noRd
apply_transformation = function(new_values, original_values) {
  unique_orig = sort(unique(original_values))
  n = length(unique_orig)
  ranks = (seq_along(unique_orig) - 1) / (n - 1)
  approx(unique_orig, ranks, xout = new_values, rule = 2)$y
}

#' Obtain an interpolation data frame
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param coords Coordinate matrix with samples as rows, and columns `x` and `y`. Rownames of coords should match colnames of mtx.
#' @param zone_obj Calibrated Zonation Object
#' @param resolution Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
#' @return An interpolation data frame
#' @noRd
apply_interpolation = function(mtx, coords, zone_obj, resolution = 1) {
  coords$zonation = getZonationGradient(mtx, zone_obj)
  nx = abs((range(coords$x)[[1]] - range(coords$x)[[2]]) / 300 * resolution)
  ny = abs((range(coords$y)[[1]] - range(coords$y)[[2]]) / 300 * resolution)
  with(coords, akima::interp(x, y, zonation, duplicate = "mean", linear = TRUE, extrap = FALSE, nx = nx, ny = ny))
}

#' Calibrate the model to baseline liver zonation
#'
#' @param mtx Gene expression matrix with genes as rows
#' coords (Optional) For spatial data, calibrates to the dimensions of the baseline sample. Coordinate matrix with samples as rows, and columns `x` and `y`.
#' @param species (Optional) Species to use, defaults to human. Currently supports 'mouse' and 'human'.
#' @return A \code{ZonationObject} with calibrated baseline zonation
#' @export
setBaseline = function(mtx, coords = NULL, species = 'human') {
  if (species == 'human') {
    hep_zonated = data.frame(read.csv(system.file('extdata', 'visium_human_regions.csv', package = 'lobular'), row.names = 1))
  } else if (species == 'mouse') {
    hep_zonated = data.frame(read.csv(system.file('extdata', 'visium_mouse_regions.csv', package = 'lobular'), row.names = 1))
  } else {
    stop("Only 'human' and 'mouse' species are supported at the moment. (Specify with species = 'mouse'")
  }

  hep_zonated = hep_zonated[rownames(hep_zonated) %in% rownames(mtx),]
  factors.zonated = hep_zonated$zonation
  names(factors.zonated) = rownames(hep_zonated)
  factors.zonated.strong = factors.zonated

  scale_factor = 1
  if (!is.null(coords)) {
    x_range = abs(range(coords$x)[[1]] - range(coords$x)[[2]])
    if (x_range < 100) {
      scale_factor = 50
    }
  }

  facts.pv = abs(factors.zonated.strong[factors.zonated.strong < 0])
  facts.cv = abs(factors.zonated.strong[factors.zonated.strong > 0])
  zonescore_raw = getGeneAvg(mtx, factors.zonated)
  zonescore_raw.pv = getGeneAvg(mtx, facts.pv)
  zonescore_raw.cv = getGeneAvg(mtx, facts.cv)
  zone_obj = ZonationObject(
    baseline = zonescore_raw,
    baseline_pv = zonescore_raw.pv,
    baseline_cv = zonescore_raw.cv,
    species = species,
    factors = factors.zonated,
    factors_pv = facts.pv,
    factors_cv = facts.cv,
    scale_factor = scale_factor
  )
  zone_obj
}

#' Apply the model to new samples, returning the zonation per cell/spot as discrete bins (1, 2, or 3)
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A vector of zonation assignments (discrete)
#' @export
getZone = function(mtx, zone_obj) {
  new_vec.pv = getGeneAvg(mtx, zone_obj$factors_pv)
  new_vec.cv = getGeneAvg(mtx, zone_obj$factors_cv)
  zone_continuous.pv = apply_transformation(new_vec.pv, zone_obj$baseline_pv)
  zone_continuous.cv = apply_transformation(new_vec.cv, zone_obj$baseline_cv)
  zonescore = (zone_continuous.cv + (1 - zone_continuous.pv)) / 2
  zone = ifelse(zonescore < (1 / 3), 'Zone_1', ifelse(zonescore < (2 / 3), 'Zone_2', 'Zone_3'))
  zone = factor(zone, levels = ZONES)
  names(zone) = colnames(mtx)
  zone
}

#' Apply the model to new samples, returning the zonation per cell/spot as a gradient between 1 and 3
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A vector of numeric zonation assignments (continuous)
#' @export
getZonationGradient = function(mtx, zone_obj) {
  new_vec = getGeneAvg(mtx, zone_obj$factors)
  zone_continuous = apply_transformation(new_vec, zone_obj$baseline)
  names(zone_continuous) = colnames(mtx)
  zone_continuous * 2 + 1
}

#' Apply the model to new samples, returning the 2d zonation per cell/spot
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A dataframe with columns PV, CV, and zone
#' @export
getZonation2d = function(mtx, zone_obj) {
  new_vec.pv = getGeneAvg(mtx, zone_obj$factors_pv)
  new_vec.cv = getGeneAvg(mtx, zone_obj$factors_cv)
  zone_continuous.pv = apply_transformation(new_vec.pv, zone_obj$baseline_pv)
  zone_continuous.cv = apply_transformation(new_vec.cv, zone_obj$baseline_cv)
  zonescore = (zone_continuous.cv + (1 - zone_continuous.pv)) / 2
  zone = ifelse(zonescore < (1 / 3), 'Zone_1', ifelse(zonescore < (2 / 3), 'Zone_2', 'Zone_3'))
  zone = factor(zone, levels = ZONES)
  zone.2d = data.frame('PV' = zone_continuous.pv, 'CV' = zone_continuous.cv, 'zone' = zone)
  rownames(zone.2d) = colnames(mtx)
  zone.2d
}

#' Apply the model to new samples, returning a plot of the 2d zonation per cell/spot
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @return A ggplot object
#' @export
plotZonation2d = function(mtx, zone_obj, point_size = 1) {
  zone.2d = getZonation2d(mtx, zone_obj)
  ggplot(zone.2d) + geom_point(aes(PV, CV, color = zone), size = point_size) + coord_fixed() +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
    geom_abline(slope = 1, intercept = 1/3, color = '#D72000FF') +
    geom_abline(slope = 1, intercept = -1/3, color = '#1BB6AFFF') +
    scale_color_manual(values = c('#1BB6AFFF', '#FFAD0AFF', '#D72000FF')) +
    labs(x = 'Portal Score', y = 'Central Score', color = 'Zone') +
    theme_classic()
}

#' Apply the model to new samples, returning a density plot of the 2d zonation per cell/spot. For visualizing zonation polarity.
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A ggplot object
#' @export
plotPolarity = function(mtx, zone_obj) {
  zone.2d = getZonation2d(mtx, zone_obj)
  polarity = -1 * cor(zone.2d$PV, zone.2d$CV)
  bins = 9
  brks = seq(0, 1, length.out = bins + 1)
  xmids = head(brks, -1) + diff(brks)/2
  ymids = xmids
  grid_counts = zone.2d %>%
    mutate(
      xbin = cut(PV, brks, include.lowest = TRUE, right = FALSE, labels = FALSE),
      ybin = cut(CV, brks, include.lowest = TRUE, right = FALSE, labels = FALSE),
      xbin = factor(xbin, levels = 1:bins),
      ybin = factor(ybin, levels = 1:bins)
    ) %>%
    count(xbin, ybin, .drop = FALSE) %>%
    mutate(
      x = xmids[as.integer(xbin)],
      y = ymids[as.integer(ybin)]
    )
  ggplot(grid_counts, aes(x, y, fill = n)) +
    geom_tile(width = 1/bins, height = 1/bins) +
    scale_fill_viridis_c(option = 'magma') +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(x = 'Portal Score', y = 'Central Score', fill = 'Count') +
    geom_abline(slope = -1 * polarity, intercept = 0.5 + polarity * 0.5, color = '#fcfdbf') +
    geom_text(aes(x = 0.5, y = 0.5, label = round(polarity, 2)),
                angle = atan(-1 * polarity) * 180 / pi, vjust = -1, color = '#fcfdbf', size = 8) +
    theme_minimal()
}

#' Use this function to infer zonation of all cell types based on their spatial interpolation within the zonation gradient.
#' Zonation gradient is inferred from a subset of samples (if specified; e.g. hepatocytes) and then applied to all samples via interpolation.
#' This function will return an output very similar to getZone, but may be slightly "smoothened".
#' An additional benefit is that this can be used to infer zonation from hepatocytes, then applied to non-parenchymal cells based on proximity to zonated hepatocytes, rather than their own gene expresison.
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param coords Coordinate matrix with samples as rows, and columns `x` and `y`. Rownames of coords should match colnames of mtx.
#' @param zone_obj Calibrated Zonation Object
#' @param resolution Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
#' @param use_for_inference (optional) A vector of sample names which should be used for zonation inference (recommended to use only hepatocytes, if annotation is available). If not provided, all samples will be used.
#' @return A vector of zonation assignments (discrete) for all samples
#' @export
getZoneSpatial = function(mtx, coords, zone_obj, resolution = 1, use_for_inference = NULL) {
  coords = data.frame(coords)
  scale_factor = zone_obj$scale_factor
  if (scale_factor > 1) {
    # Move into a reasonable range for interpolation
    coords$x = coords$x * scale_factor
    coords$y = coords$y * scale_factor
    resolution = resolution * scale_factor
  }
  if (length(use_for_inference)) {
    coords_subset = coords[rownames(coords) %in% use_for_inference,]
    mtx = mtx[,colnames(mtx) %in% use_for_inference]
  } else {
    coords_subset = coords
  }
  interp_data = apply_interpolation(mtx, coords_subset, zone_obj, resolution)
  ix = findInterval(coords$x, interp_data$x)
  iy = findInterval(coords$y, interp_data$y)
  ix[ix == 0] = 1
  iy[iy == 0] = 1
  ix[ix == length(interp_data$x)] = length(interp_data$x) - 1
  iy[iy == length(interp_data$y)] = length(interp_data$y) - 1
  interp_value = mapply(function(i, j, x0, y0) {
    x1 = interp_data$x[i]; x2 = interp_data$x[i + 1]
    y1 = interp_data$y[j]; y2 = interp_data$y[j + 1]
    z11 = interp_data$z[i, j]; z21 = interp_data$z[i + 1, j]
    z12 = interp_data$z[i, j + 1]; z22 = interp_data$z[i + 1, j + 1]
    if (any(is.na(c(z11, z21, z12, z22)))) return(NA)
    wx = (x0 - x1) / (x2 - x1)
    wy = (y0 - y1) / (y2 - y1)
    z = (1 - wx) * (1 - wy) * z11 + wx * (1 - wy) * z21 +
      (1 - wx) * wy * z12 + wx * wy * z22
    return(z)
  }, ix, iy, coords$x, coords$y)
  coords$zonation = interp_value
  breaks = seq(1, 3, length.out = 4)
  zone = cut(coords$zonation, breaks = breaks, labels = ZONES, include.lowest = TRUE)
  zone = factor(zone, levels = ZONES)
  names(zone) = rownames(coords)
  zone
}

#' Plots the 2-dimensional interpolated zones in a spatial dataset
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param coords Coordinate matrix with samples as rows, and columns `x` and `y`. Rownames of coords should match colnames of mtx.
#' @param zone_obj Calibrated Zonation Object
#' @param resolution Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
#' @param use_for_inference (optional) A vector of sample names which should be used for zonation inference (recommended to use only hepatocytes, if annotation is available). If not provided, all samples will be used.
#' @return A ggplot object
#' @export
plotZoneSpatial = function(mtx, coords, zone_obj, resolution = 1, use_for_inference = NULL) {
  coords = data.frame(coords)
  scale_factor = zone_obj$scale_factor
  if (scale_factor > 1) {
    # Move into a reasonable range for interpolation
    coords$x = coords$x * scale_factor
    coords$y = coords$y * scale_factor
    resolution = resolution * scale_factor
  }
  if (length(use_for_inference)) {
    coords_subset = coords[rownames(coords) %in% use_for_inference,]
    mtx = mtx[,colnames(mtx) %in% use_for_inference]
  } else {
    coords_subset = coords
  }
  interp_data = apply_interpolation(mtx, coords_subset, zone_obj, resolution)
  interp_df = data.frame(
    x = rep(interp_data$x, times = length(interp_data$y)),
    y = rep(interp_data$y, each = length(interp_data$x)),
    z = as.vector(interp_data$z)
  )
  if (scale_factor > 1) {
    # Transfer back to original space
    interp_df$x = interp_df$x / scale_factor
    interp_df$y = interp_df$y / scale_factor
  }
  breaks = seq(1, 3, length.out = 4)
  ggplot(interp_df, aes(x, y, z = z)) +
    geom_contour_filled(breaks = breaks) +
    coord_fixed() +
    scale_fill_manual(name = 'Zone', labels = paste('Zone', 1:3), values = c('#333F48', '#C6AA76', '#BA0C2F')) # "Bull City" color scheme from https://r-graph-gallery.com/color-palette-finder
}

#' Plots the 2-dimensional interpolated zones with zonation contour outlines in a spatial dataset
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param coords Coordinate matrix with samples as rows, and columns `x` and `y`. Rownames of coords should match colnames of mtx.
#' @param zone_obj Calibrated Zonation Object
#' @param resolution Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @param use_for_inference (optional) A vector of sample names which should be used for zonation inference (recommended to use only hepatocytes, if annotation is available). If not provided, all samples will be used.
#' @return A ggplot object
#' @export
plotZoneSpatialContours = function(mtx, coords, zone_obj, resolution = 1, point_size = 1, plot_options = NULL, use_for_inference = NULL) {
  coords = data.frame(coords)
  scale_factor = zone_obj$scale_factor
  if (scale_factor > 1) {
    # Move into a reasonable range for interpolation
    coords$x = coords$x * scale_factor
    coords$y = coords$y * scale_factor
    resolution = resolution * scale_factor
  }
  if (length(use_for_inference)) {
    coords_subset = coords[rownames(coords) %in% use_for_inference,]
    mtx_subset = mtx[,colnames(mtx) %in% use_for_inference]
  } else {
    coords_subset = coords
    mtx_subset = mtx
  }
  interp_data = apply_interpolation(mtx_subset, coords_subset, zone_obj, resolution)
  zone_assignments = getZonationGradient(mtx, zone_obj)
  coords$zone = zone_assignments
  interp_df = data.frame(
    x = rep(interp_data$x, times = length(interp_data$y)),
    y = rep(interp_data$y, each = length(interp_data$x)),
    z = as.vector(interp_data$z)
  )
  if (scale_factor > 1) {
    # Transfer back to original space
    coords$x = coords$x / scale_factor
    coords$y = coords$y / scale_factor
    interp_df$x = interp_df$x / scale_factor
    interp_df$y = interp_df$y / scale_factor
  }
  breaks = seq(1, 3, length.out = 4)
  ggplot(coords) +
    geom_point(data = coords, aes(x = x, y = y, color = zone), size = point_size) +
    scale_color_viridis_c(name = 'Zonation') +
    geom_contour(data = interp_df,
                aes(x = x, y = y, z = z),
                breaks = breaks,
                color = 'black',
                linewidth = 2) +
    coord_fixed()
}

#' Plots a custom variable with zonation contour outlines in a spatial dataset
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param meta Metadata matrix with samples as rows, and columns `x`, `y`, and `mycolname`, where `mycolname` is passed as `colname`. Rownames of coords should match colnames of mtx.
#' @param colname Name of custom column in `meta`
#' @param zone_obj Calibrated Zonation Object
#' @param resolution Optional numeric value for the resolution, where higher value results in a more granular interpolation (default 1)
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @param use_for_inference (optional) A vector of sample names which should be used for zonation inference (recommended to use only hepatocytes, if annotation is available). If not provided, all samples will be used.
#' @return A ggplot object
#' @export
plotZoneSpatialCustom = function(mtx, meta, col_name, zone_obj, resolution = 1, point_size = 1, use_for_inference = NULL) {
  meta = data.frame(meta)
  scale_factor = zone_obj$scale_factor
  if (scale_factor > 1) {
    # Move into a reasonable range for interpolation
    meta$x = meta$x * scale_factor
    meta$y = meta$y * scale_factor
    resolution = resolution * scale_factor
  }
  if (length(use_for_inference)) {
    meta_subset = meta[rownames(meta) %in% use_for_inference,]
    mtx_subset = mtx[,colnames(mtx) %in% use_for_inference]
  } else {
    meta_subset = meta
    mtx_subset = mtx
  }
  interp_data = apply_interpolation(mtx_subset, meta_subset, zone_obj, resolution)
  interp_df = data.frame(
    x = rep(interp_data$x, times = length(interp_data$y)),
    y = rep(interp_data$y, each = length(interp_data$x)),
    z = as.vector(interp_data$z)
  )
  if (scale_factor > 1) {
    # Transfer back to original space
    meta$x = meta$x / scale_factor
    meta$y = meta$y / scale_factor
    interp_df$x = interp_df$x / scale_factor
    interp_df$y = interp_df$y / scale_factor
  }
  breaks = seq(1, 3, length.out = 4)
  ggplot(meta) +
    geom_point(data = meta, aes(x = x, y = y, color = .data[[col_name]]), size = point_size) +
    geom_contour(data = interp_df,
                aes(x = x, y = y, z = z),
                breaks = breaks,
                color = 'black',
                linewidth = 2) +
    coord_fixed()
}
