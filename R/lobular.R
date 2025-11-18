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
ZonationObject = function(baseline_1, baseline_2, baseline_3, species, factors_1, factors_2, factors_3, scale_factor) {
  # Input validation
  if (!species %in% c("human", "mouse")) {
    stop("species must be 'human' or 'mouse'")
  }
  if (!is.numeric(scale_factor) || is.null(scale_factor)) {
    stop("scale_factor must be a number")
  }

  # Create the object
  obj = list(
    baseline_1 = baseline_1,
    baseline_2 = baseline_2,
    baseline_3 = baseline_3,
    species = species,
    factors_1 = factors_1,
    factors_2 = factors_2,
    factors_3 = factors_3,
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
#' @param factor_threshold (Optional) Minimum value for zonation factors to be included in calculation (removes noise).
#' @return A \code{ZonationObject} with calibrated baseline zonation
#' @export
setBaseline = function(mtx, coords = NULL, species = 'human', factor_threshold = 0.1) {
  if (species == 'human') {
    orig_annotation = hep_zonated = data.frame(read.csv(system.file('extdata', 'visium_human_zonation.csv', package = 'lobular'), row.names = 1))
  } else if (species == 'mouse') {
    orig_annotation = hep_zonated = data.frame(read.csv(system.file('extdata', 'visium_mouse_zonation.csv', package = 'lobular'), row.names = 1))
  } else {
    stop("Only 'human' and 'mouse' species are supported at the moment. (Specify with species = 'mouse'")
  }

  minMaxNorm = function(v) (v - min(v, na.rm = T)) / (max(v, na.rm = T) - min(v, na.rm = T))
  getZoneObj = function(df) {
    factors.zonated.1 = df$zone_1
    factors.zonated.2 = df$zone_2
    factors.zonated.3 = df$zone_3
    names(factors.zonated.1) = rownames(df)
    names(factors.zonated.2) = rownames(df)
    names(factors.zonated.3) = rownames(df)
    factors.zonated.1 = ifelse(is.na(factors.zonated.1), 0, factors.zonated.1^2)
    factors.zonated.2 = ifelse(is.na(factors.zonated.2), 0, factors.zonated.2^2)
    factors.zonated.3 = ifelse(is.na(factors.zonated.3), 0, factors.zonated.3^2)
    factors.zonated.1 = ifelse(factors.zonated.1 > factor_threshold, factors.zonated.1, 0)
    factors.zonated.2 = ifelse(factors.zonated.2 > factor_threshold, factors.zonated.2, 0)
    factors.zonated.3 = ifelse(factors.zonated.3 > factor_threshold, factors.zonated.3, 0)

    scale_factor = 1
    if (!is.null(coords)) {
      x_range = abs(range(coords$x)[[1]] - range(coords$x)[[2]])
      if (x_range < 100) {
        scale_factor = 50
      }
    }

    zonescore_raw.1 = getGeneAvg(mtx, factors.zonated.1)
    zonescore_raw.2 = getGeneAvg(mtx, factors.zonated.2)
    zonescore_raw.3 = getGeneAvg(mtx, factors.zonated.3)

    ZonationObject(
      baseline_1 = zonescore_raw.1,
      baseline_2 = zonescore_raw.2,
      baseline_3 = zonescore_raw.3,
      species = species,
      factors_1 = factors.zonated.1,
      factors_2 = factors.zonated.2,
      factors_3 = factors.zonated.3,
      scale_factor = scale_factor
    )
  }

  hep_zonated = hep_zonated[rownames(hep_zonated) %in% rownames(mtx),]
  zone_obj = getZoneObj(hep_zonated)

  grad = getZonationGradient(mtx, zone_obj)
  grad = grad[colnames(mtx)]
  zone2_grad = minMaxNorm(grad)
  zone2_grad = ifelse(zone2_grad < 0.5, zone2_grad * 2, (1 - zone2_grad) * 2)
  cor_zone2_grad = apply(mtx, 1, function(row) cor(zone2_grad, row))
  cor_grad = apply(mtx, 1, function(row) cor(grad, row))
  cor_grad_1 = ifelse(cor_grad < 0, abs(cor_grad), 0)
  cor_grad_2 = ifelse(cor_zone2_grad > 0, cor_zone2_grad, 0)
  cor_grad_3 = ifelse(cor_grad > 0, cor_grad, 0)
  hep_zonated = data.frame(zone_1 = cor_grad_1, zone_2 = cor_grad_2, zone_3 = cor_grad_3)
  zone_obj_2 = getZoneObj(hep_zonated)

  getSimilarity = function(orig, curr, zone) {
    orig = orig[order(orig[[zone]], decreasing = T),]^2
    orig_genes = rownames(orig)[1:10]
    curr_vals = curr[orig_genes, zone]
    curr_vals = ifelse(is.na(curr_vals), 0, curr_vals)
    print(paste0('This baseline is ', round(cor(orig[orig_genes, zone], curr_vals, method = 'pearson'), 2) * 100, '% similar to the reference (based on top 10 reference genes).'))
  }
  getSimilarity(orig_annotation, hep_zonated, 'zone_1')

  zone_obj_2
}

#' Get the pearson correlations between zone scores and genes
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A data frame containing the pearson correlations between zone scores and genes
#' @export
getGeneZonation = function(zone_obj) {
  sqrt(data.frame(zone_1 = zone_obj$factors_1, zone_2 = zone_obj$factors_2, zone_3 = zone_obj$factors_3))
}

#' Apply the model to new samples, returning the zonation per cell/spot as a gradient between 1 and 3
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A vector of numeric zonation assignments (continuous)
#' @export
getZonationGradient = function(mtx, zone_obj) {
  new_vec.1 = getGeneAvg(mtx, zone_obj$factors_1)
  new_vec.3 = getGeneAvg(mtx, zone_obj$factors_3)
  zone_continuous.1 = apply_transformation(new_vec.1, zone_obj$baseline_1)
  zone_continuous.3 = apply_transformation(new_vec.3, zone_obj$baseline_3)
  zonescore = (zone_continuous.3 + (1 - zone_continuous.1)) / 2
  names(zonescore) = colnames(mtx)
  zonescore * 2 + 1
}

#' Apply the model to new samples, returning the zonation per cell/spot as discrete bins (1, 2, or 3)
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A vector of zonation assignments (discrete)
#' @export
getZone = function(mtx, zone_obj) {
  zonescore = getZonationGradient(mtx, zone_obj)
  zone = ifelse(zonescore < 1 * (2 / 3) + 1, 'Zone_1', ifelse(zonescore < 2 * (2 / 3) + 1, 'Zone_2', 'Zone_3'))
  zone = factor(zone, levels = ZONES)
  names(zone) = colnames(mtx)
  zone
}

#' Apply the model to new samples, returning the 2d zonation per cell/spot
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A dataframe with columns PV, CV, and zone
#' @export
getZonation2d = function(mtx, zone_obj) {
  new_vec.1 = getGeneAvg(mtx, zone_obj$factors_1)
  new_vec.3 = getGeneAvg(mtx, zone_obj$factors_3)
  zone_continuous.1 = apply_transformation(new_vec.1, zone_obj$baseline_1)
  zone_continuous.3 = apply_transformation(new_vec.3, zone_obj$baseline_3)
  if (zone_obj$species == 'mouse') {
    new_vec.2 = getGeneAvg(mtx, zone_obj$factors_2)
    zone_continuous.2 = apply_transformation(new_vec.2, zone_obj$baseline_2)
  } else {
    zone_continuous.2 = rep(0, length(zone_continuous.1))
  }
  zonescore = (zone_continuous.3 + (1 - zone_continuous.1)) / 2
  zone = ifelse(zonescore < (1 / 3), 'Zone_1', ifelse(zonescore < (2 / 3), 'Zone_2', 'Zone_3'))
  zone = factor(zone, levels = ZONES)
  zone.2d = data.frame('ZONE_1' = zone_continuous.1, 'ZONE_2' = zone_continuous.2, 'ZONE_3' = zone_continuous.3, 'zone' = zone)
  rownames(zone.2d) = colnames(mtx)
  zone.2d
}

#' Apply the model to new samples, returning a plot of the 2d zonation per cell/spot with *zone* indicated by the color
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @return A ggplot object
#' @export
plotZonation2d = function(mtx, zone_obj, point_size = 1) {
  zone.2d = getZonation2d(mtx, zone_obj)
  ggplot(zone.2d) + geom_point(aes(ZONE_1, ZONE_3, color = zone), size = point_size) + coord_fixed() +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
    geom_abline(slope = 1, intercept = 1/3, color = '#D72000FF') +
    geom_abline(slope = 1, intercept = -1/3, color = '#1BB6AFFF') +
    scale_color_manual(values = c('#1BB6AFFF', '#FFAD0AFF', '#D72000FF')) +
    labs(x = 'Zone 1 Score', y = 'Zone 3 Score', color = 'Zone') +
    ggdark::dark_theme_classic()
}

#' Apply the model to new samples, returning a plot of the 2d zonation per cell/spot with *zone 2 score* indicated by the color. Only supported for mouse (no human zone 2 genes).
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @return A ggplot object
#' @export
plotZonation2d_2 = function(mtx, zone_obj, point_size = 1) {
  if (zone_obj$species != 'mouse') {
    stop(paste0("This plot is not available for ", zone_obj$species, ", as there are no zone 2 reference genes."))
  }
  zone.2d = getZonation2d(mtx, zone_obj)
  ggplot(zone.2d) + geom_point(aes(ZONE_1, ZONE_3, color = ZONE_2), size = point_size) + coord_fixed() +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
    geom_abline(slope = 1, intercept = 1/3, color = '#FFFFFF') +
    geom_abline(slope = 1, intercept = -1/3, color = '#FFFFFF') +
    scale_color_viridis_c() +
    labs(x = 'Zone 1 Score', y = 'Zone 3 Score', color = 'Zone 2 Score') +
    ggdark::dark_theme_classic()
}

#' Apply the model to new samples, returning a ridge plot of the scores for zones 1, 2, and 3
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A ggplot object
#' @export
plotZonationRidge = function(mtx, zone_obj) {
  zone.2d = getZonation2d(mtx, zone_obj)
  zone.2d.long <- zone.2d[,1:3] %>%
    tidyr::pivot_longer(cols = c(ZONE_1, ZONE_2, ZONE_3),
                  names_to = 'Zone',
                  values_to = 'Value')
  zone.2d.long$Zone = factor(zone.2d.long$Zone, levels = c('ZONE_3', 'ZONE_2', 'ZONE_1'))
  ggplot(zone.2d.long, aes(x = Value, y = Zone, fill = Zone)) +
    ggridges::geom_density_ridges(scale = 1.2, rel_min_height = 0.01) +
    scale_fill_manual(values = c('#1BB6AFFF', '#FFAD0AFF', '#D72000FF')) +
    ggdark::dark_theme_classic() +
    theme(legend.position = 'none')
}

#' Apply the model to new samples, returning a density plot of the 2d zonation per cell/spot. For visualizing zonation polarity.
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A ggplot object
#' @export
plotPolarity = function(mtx, zone_obj) {
  zone.2d = getZonation2d(mtx, zone_obj)
  polarity = -1 * cor(zone.2d$ZONE_1, zone.2d$ZONE_3)
  bins = 9
  brks = seq(0, 1, length.out = bins + 1)
  xmids = head(brks, -1) + diff(brks)/2
  ymids = xmids
  grid_counts = zone.2d %>%
    dplyr::mutate(
      xbin = cut(ZONE_1, brks, include.lowest = TRUE, right = FALSE, labels = FALSE),
      ybin = cut(ZONE_3, brks, include.lowest = TRUE, right = FALSE, labels = FALSE),
      xbin = factor(xbin, levels = 1:bins),
      ybin = factor(ybin, levels = 1:bins)
    ) %>%
    dplyr::count(xbin, ybin, .drop = FALSE) %>%
    dplyr::mutate(
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
    ggdark::dark_theme_classic()
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
#' @param line_width Optional numeric value for the ggplot contour line width (default 2)
#' @param use_for_inference (optional) A vector of sample names which should be used for zonation inference (recommended to use only hepatocytes, if annotation is available). If not provided, all samples will be used.
#' @return A ggplot object
#' @export
plotZoneSpatialContours = function(mtx, coords, zone_obj, resolution = 1, point_size = 1, line_width = 2, plot_options = NULL, use_for_inference = NULL) {
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
                linewidth = line_width) +
    coord_fixed() +
    dark_theme_classic()
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
