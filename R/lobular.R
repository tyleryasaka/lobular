ZONES = c('Zone_1', 'Zone_2', 'Zone_3')
ZONE_COLORS = c('#504880FF', '#1BB6AFFF', '#FFAD0AFF', '#D72000FF', '#F080D8FF')
minMaxNorm = function(v) (v - min(v, na.rm = T)) / (max(v, na.rm = T) - min(v, na.rm = T))
.lobular_cache = new.env(parent = emptyenv())
.obj_id_counter = new.env(parent = emptyenv())
.obj_id_counter$n = 0L

next_obj_id = function() {
  .obj_id_counter$n <- .obj_id_counter$n + 1L
  as.character(.obj_id_counter$n)
}

access_cache = function(zone_obj, prop, fn) {
  id = zone_obj$obj_id
  if (is.null(.lobular_cache[[id]])) {
    .lobular_cache[[id]] = new.env(parent = emptyenv())
  }
  if (is.null(.lobular_cache[[id]][[prop]])) {
    .lobular_cache[[id]][[prop]] = fn()
  }
  .lobular_cache[[id]][[prop]]
}

normalizeMatrix = function(mtx) {
  mtx = as.matrix(mtx)
  mtx[is.na(mtx)] = 0
  mtx
}

getZonationCor = function(mtx, zonation) {
  cor_results <- t(apply(mtx, 1, function(x) {
    test <- cor.test(x, zonation, method = "spearman", exact = FALSE)
    c(rho = test$estimate, pval = test$p.value)
  }))
  cor_results <- as.data.frame(cor_results)
  cor_results$padj <- p.adjust(cor_results$pval, method = "BH")
  cor_results <- cor_results[order(cor_results$padj), ]
  cor_results
}

check_input_consistency = function(model) {
  mtx = model$mtx
  common = intersect(rownames(mtx), names(model$gene_means))
  if (length(common) == 0) stop("No genes in common between mtx and model")
  input_means = rowMeans(mtx[common, , drop = FALSE])
  stored_means = model$gene_means[common]
  ratio = median(abs(input_means), na.rm = TRUE) / median(abs(stored_means), na.rm = TRUE)
  if (ratio > 10 || ratio < 0.1)
    warning(sprintf("Input gene means differ from training gene means by ~%.0fx. Check that mtx is on the expected scale.", ratio))
}

predict_position = function(model) {
  access_cache(model, 'predict_position', function() {
    mtx = model$mtx
    check_input_consistency(model)
    pred_mtx = mtx
    common_genes = intersect(rownames(pred_mtx), names(model$weights))
    pred_mtx = pred_mtx[common_genes, , drop = FALSE]
    pred_mtx = scale(t(pred_mtx), center = model$gene_means[common_genes], scale = FALSE)
    w_s = model$weights[common_genes]
    if(length(w_s) > 0) w_s = w_s / sqrt(sum(w_s^2))
    raw = as.vector(pred_mtx %*% w_s)
    (model$f_e(raw) - model$f_rl) / (model$f_rh - model$f_rl)
  })
}

predict_position_2d = function(model) {
  access_cache(model, 'predict_position_2d', function() {
    mtx = model$mtx
    pg = predict_position(model)
    w_sub = if (!is.null(model$weights_unfilt)) model$weights_unfilt else model$weights
    common = intersect(rownames(mtx), names(w_sub))
    mtx_c = scale(t(mtx[common, , drop = FALSE]), center = model$gene_means[common], scale = FALSE)
    gp = function(sub_m) {
      valid_gs = intersect(sub_m$genes, common)
      ws = w_sub[common]
      ws[!(names(ws) %in% valid_gs)] = 0
      raw = as.vector(mtx_c %*% ws)
      vals = sub_m$e(raw)
      (vals - sub_m$e_min) / (sub_m$e_max - sub_m$e_min)
    }
    z = ifelse(pg > model$q2, 'Zone_3', ifelse(pg > model$q1, 'Zone_2', 'Zone_1'))
    data.frame(ZONE_1 = 1 - gp(model$m1), ZONE_3 = gp(model$m3),
               zonation = pg, zone = factor(z, levels = c('Zone_1', 'Zone_2', 'Zone_3')))
  })
}

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
apply_transformation = function(new_values, original_values, p = 0) {
  lims = quantile(original_values, c(p, 1 - p))
  orig = pmax(pmin(original_values, lims[2]), lims[1])
  u = sort(unique(orig))
  n = length(u)
  approx(u, (seq_along(u) - 1) / max(1, n - 1), xout = new_values, rule = 2)$y
}

getZonationGradient_help = function(zone_obj) {
  predict_position(zone_obj) * 2 + 1
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
  coords$zonation = getZonationGradient_help(mtx, zone_obj)
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
trainModel = function(mtx, coords = NULL, species = 'human', regularization = 0.8, filter = 0, verbose = FALSE) {
  if (species == 'human') {
    initial_weights = readRDS(system.file('extdata', 'initial_weights_human.RDS', package = 'lobular'))
  } else if (species == 'mouse') {
    initial_weights = readRDS(system.file('extdata', 'initial_weights_mouse.RDS', package = 'lobular'))
  } else {
    stop("Only 'human' and 'mouse' species are supported at the moment. (Specify with species = 'mouse'")
  }
  initial_weights = head(initial_weights[order(abs(initial_weights), decreasing = T)], 500)
  mtx = normalizeMatrix(mtx)
  em_zonation(mtx, initial_weights, iterations = 10, density_cut = 0, regularization = regularization, cor_thresh = filter, verbose = verbose)
}

applyModel = function(mtx, zone_obj) {
  zone_obj$mtx = normalizeMatrix(mtx)
  zone_obj$obj_id = next_obj_id()
  zone_obj
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
getZonationGradient = function(zone_obj) {
  getZonationGradient_help(zone_obj)
}

#' Apply the model to new samples, returning the zonation per cell/spot as discrete bins (1, 2, or 3)
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A vector of zonation assignments (discrete)
#' @export
getZone = function(zone_obj) {
  mtx = zone_obj$mtx
  zonescore = predict_position(zone_obj)
  zone = ifelse(zonescore < (1 / 3), 'Zone_1', ifelse(zonescore < (2 / 3), 'Zone_2', 'Zone_3'))
  zone = factor(zone, levels = ZONES)
  names(zone) = colnames(mtx)
  zone
}

getZonation2d_help = function(zone_obj) {
  mtx = zone_obj$mtx
  zone.2d = predict_position_2d(zone_obj)
  rownames(zone.2d) = colnames(mtx)
  zone.2d
}

#' Apply the model to new samples, returning the 2d zonation per cell/spot
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A dataframe with columns PV, CV, and zone
#' @export
getZonation2d = function(zone_obj) {
  getZonation2d_help(zone_obj)
}

#' Apply the model to new samples, returning a plot of the 2d zonation per cell/spot with *zone* indicated by the color
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @return A ggplot object
#' @export
plotZonation2d = function(zone_obj, point_size = 1) {
  zone.2d = getZonation2d_help(zone_obj)
  ggplot(zone.2d) + geom_point(aes(ZONE_1, ZONE_3, color = zone), size = point_size) + coord_fixed() +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_color_manual(values = ZONE_COLORS, limits = ZONES) +
    labs(x = 'Zone 1 Score', y = 'Zone 3 Score', color = 'Zone')
}

#' Apply the model to new samples, returning a plot of the 2d zonation per cell/spot with *zone 2 score* indicated by the color. Only supported for mouse (no human zone 2 genes).
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @return A ggplot object
#' @export
plotZonation2d_2 = function(zone_obj, point_size = 1) {
  if (zone_obj$species != 'mouse') {
    stop(paste0("This plot is not available for ", zone_obj$species, ", as there are no zone 2 reference genes."))
  }
  zone.2d = getZonation2d_help(zone_obj)
  ggplot(zone.2d) + geom_point(aes(ZONE_1, ZONE_3, color = ZONE_2), size = point_size) + coord_fixed() +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
    geom_abline(slope = 1, intercept = 1/3, color = '#FFFFFF') +
    geom_abline(slope = 1, intercept = -1/3, color = '#FFFFFF') +
    scale_color_viridis_c() +
    labs(x = 'Zone 1 Score', y = 'Zone 3 Score', color = 'Zone 2 Score') +
    ggdark::dark_theme_classic()
}

#' Apply the model to new samples, returning a plot of the 2d zonation per cell/spot with *expression of a specific gene* indicated by the color
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @param gene Name of gene to plot
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @return A ggplot object
#' @export
plotZonation2dGene = function(zone_obj, gene, point_size = 1) {
  mtx = zone_obj$mtx
  zone.2d = getZonation2d_help(zone_obj)
  zone.2d[,gene] = mtx[gene,]
  ggplot(zone.2d) + geom_point(aes(ZONE_1, ZONE_3, color = .data[[gene]]), size = point_size) + coord_fixed() +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_color_viridis_c() +
    labs(x = 'Zone 1 Score', y = 'Zone 3 Score', color = gene)
}

#' Apply the model to new samples, returning a ridge plot of gene expression along zonation axis
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @param gene Name of gene to plot
#' @return A ggplot object
#' @export
plotRegression = function(zone_obj, gene) {
  mtx = zone_obj$mtx
  zone.2d = getZonation2d_help(zone_obj)
  zone.2d[,gene] = mtx[gene,]
  ggplot(zone.2d, aes(x = zonation, y = .data[[gene]])) +
    geom_point(color = '#504880FF') +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE, color = "#F080D8FF", fill  = "#C8C0F8FF") +
    xlab('Zonation')
}

#' Apply the model to new samples, returning a density plot of the 2d zonation per cell/spot. For visualizing zonation polarity.
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A ggplot object
#' @export
plotPolarity = function(zone_obj) {
  zone.2d = getZonation2d_help(zone_obj)
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
                angle = atan(-1 * polarity) * 180 / pi, vjust = -1, color = '#fcfdbf', size = 8)
}

#' Generate a plot of the "diff" or breakdown of differential expression of zonated genes between 2 samples
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone Numeric, either 1, 2 (mouse only), or 3
#' @param zone_obj Calibrated Zonation Object
#' @param threshold Threshold of zonation for genes to include in plot
#' @return A ggplot object
#' @export
plotZonationDiff = function(mtx_1, mtx_2, zone_obj, zone, threshold = 0.1, font_size = 9) {
  mtx_1 = normalizeMatrix(mtx_1)
  mtx_2 = normalizeMatrix(mtx_2)
  allowed_zones = c(1, 3)
  if (zone_obj$species == 'mouse') {
    allowed_zones = c(1, 2, 3)
  }
  if (!(zone %in% allowed_zones)) {
    stop(paste0('For ', zone_obj$species, ', `zone` must be one of: ', paste(allowed_zones, collapse = ',')))
  }
  factors_zone = zone_obj[[paste0('factors_', zone)]]
  factors_zone = factors_zone[factors_zone > threshold]
  common_genes = intersect(names(factors_zone), rownames(mtx_1))
  common_genes = intersect(common_genes, rownames(mtx_2))
  factors_zone = factors_zone[common_genes]
  factors_zone = factors_zone[order(factors_zone, decreasing = T)]
  mtx_filtered = mtx_2[common_genes,]
  bl_filtered = mtx_1[common_genes,]
  diff = Matrix::rowMeans(mtx_filtered) - Matrix::rowMeans(bl_filtered)
  diff = diff[names(factors_zone)]
  df = data.frame(id = seq_along(diff), name = names(diff), value = as.numeric(diff), weight = factors_zone)
  gap = 0.05
  df$xmin <- NA_real_
  df$xmax <- NA_real_
  df$xmin[1] <- 0
  df$xmax[1] <- df$xmin[1] + df$weight[1]
  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      df$xmin[i] <- df$xmax[i-1] + gap
      df$xmax[i] <- df$xmin[i] + df$weight[i]
    }
  }
  df$center <- (df$xmin + df$xmax)/2
  df$name = factor(df$name, levels = df$name)
  ggplot(df) +
    geom_rect(aes(
      xmin = xmin,
      xmax = xmax,
      ymin = pmin(0, value),
      ymax = pmax(0, value),
      fill = value > 0
    )) +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = c("TRUE" = "#4daf4a", "FALSE" = "#e41a1c")) +
    scale_x_continuous(
      breaks = df$center,
      labels = df$name
    ) +
    ggdark::dark_theme_classic() +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1, size = font_size)) +
    labs(x = 'Gene', y = "Expression Change", title = paste0('Differential Zonation - Zone ', zone))
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
getZoneSpatial = function(coords, zone_obj, resolution = 1, use_for_inference = NULL) {
  mtx = zone_obj$mtx
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
  zone = cut(coords$zonation, breaks = breaks, labels = c('Zone_1', 'Zone_2', 'Zone_3'), include.lowest = TRUE)
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
plotZoneSpatial = function(coords, zone_obj, resolution = 1, use_for_inference = NULL) {
  mtx = zone_obj$mtx
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
plotZoneSpatialContours = function(coords, zone_obj, resolution = 1, point_size = 1, line_width = 2, plot_options = NULL, use_for_inference = NULL) {
  mtx = zone_obj$mtx
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
plotZoneSpatialCustom = function(meta, col_name, zone_obj, resolution = 1, point_size = 1, use_for_inference = NULL) {
  mtx = zone_obj$mtx
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

#' Plot a virtual lobule of the inferred zonation distribution
#'
#' Visualizes the distribution of inferred zonation within an idealized
#' hexagonal lobule. Each pixel is placed according to its rank along the
#' corner-to-center axis (0 = nearest corner / Zone 1, 1 = center / Zone 3)
#' and colored by the matching quantile of \code{getZonationGradient(zone_obj)}.
#'
#' @param zone_obj A calibrated Zonation Object (output of \code{applyModel}).
#' @param resolution Integer pixel resolution along the x-axis (default 100).
#' @param palette Character; a \pkg{paletteer} continuous palette name
#' @param reverse_palette Logical; reverse palette direction (default FALSE).
#' @param pointy_top Logical; if TRUE the hexagon has a vertex at the top,
#'   otherwise it is flat-topped (default FALSE, matching \code{virtual_lobule}).
#' @param show_legend Logical; show the colour legend (default TRUE).
#' @param seed Optional integer for reproducibility. Not strictly required
#'   here (the mapping is deterministic), but kept for API symmetry with
#'   \code{virtual_lobule}.
#' @return A ggplot object.
#' @export
plotVirtualLobule <- function(zone_obj,
                              resolution = 100,
                              palette = "ggthemes::Classic Red-Blue",
                              reverse_palette = T,
                              pointy_top = FALSE,
                              show_legend = TRUE,
                              seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # --- 0. Pull the inferred zonation gradient (values on [1, 3]) ---
  zonation <- getZonationGradient(zone_obj)
  zonation <- zonation[is.finite(zonation)]
  if (length(zonation) == 0) {
    stop("No finite zonation values from getZonationGradient(zone_obj).")
  }

  # --- 1. Hexagon geometry ---
  R <- 1
  angle_offset <- if (pointy_top) pi / 6 else 0
  vertex_angles <- angle_offset + (0:5) * pi / 3
  vertices_x <- R * cos(vertex_angles)
  vertices_y <- R * sin(vertex_angles)

  # --- 2. Pixel grid ---
  x_range <- range(vertices_x)
  y_range <- range(vertices_y)
  pixel_size <- diff(x_range) / resolution

  x_seq <- seq(x_range[1] + pixel_size / 2, x_range[2] - pixel_size / 2, by = pixel_size)
  y_seq <- seq(y_range[1] + pixel_size / 2, y_range[2] - pixel_size / 2, by = pixel_size)
  grid <- expand.grid(x = x_seq, y = y_seq)

  # --- 3. Point-in-hexagon test (cross-product against CCW edges) ---
  inside <- rep(TRUE, nrow(grid))
  for (i in 1:6) {
    j <- (i %% 6) + 1
    ex <- vertices_x[j] - vertices_x[i]
    ey <- vertices_y[j] - vertices_y[i]
    cross <- ex * (grid$y - vertices_y[i]) - ey * (grid$x - vertices_x[i])
    inside <- inside & (cross >= 0)
  }
  grid <- grid[inside, , drop = FALSE]
  if (nrow(grid) == 0) stop("No pixels inside hexagon. Increase resolution.")

  # --- 4. Position: 0 = nearest corner (portal), 1 = center (central) ---
  d_center <- sqrt(grid$x^2 + grid$y^2)
  d_nearest_corner <- do.call(pmin, lapply(seq_len(6), function(i) {
    sqrt((grid$x - vertices_x[i])^2 + (grid$y - vertices_y[i])^2)
  }))
  raw_pos <- d_nearest_corner / (d_nearest_corner + d_center)
  raw_pos[d_center < 1e-12] <- 1

  # Rank-normalise pixel positions to a uniform [0, 1] distribution
  r <- rank(raw_pos, ties.method = "average")
  grid$position <- (r - min(r)) / (max(r) - min(r))

  # --- 5. Map each pixel to the matching quantile of the zonation gradient ---
  grid$zonation <- as.numeric(stats::quantile(
    zonation,
    probs = grid$position,
    na.rm = TRUE,
    names = FALSE,
    type  = 7
  ))

  # --- 6. Colour palette ---
  n_colors <- 256
  pal_colors <- tryCatch(
    as.character(paletteer::paletteer_c(palette, n = n_colors)),
    error = function(e) as.character(paletteer::paletteer_d(palette))
  )
  if (reverse_palette) pal_colors <- rev(pal_colors)

  # --- 7. ggplot ---
  ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, fill = zonation)) +
    ggplot2::geom_tile(width = pixel_size, height = pixel_size) +
    ggplot2::scale_fill_gradientn(
      colours  = pal_colors,
      limits   = c(1, 3),
      oob      = scales::squish,
      na.value = "grey50",
      name     = "Zonation"
    ) +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position  = if (show_legend) "right" else "none",
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA)
    )
}

#' Density plot of inferred zonation across one or more samples
#'
#' Overlays the distributions of \code{getZonationGradient()} values for one
#' or more samples as outlined density curves on a shared baseline.
#'
#' @param zone_objs A single Zonation Object, or a (preferably named) list of
#'   them. List names are used as legend labels.
#' @param palette Character; a paletteer palette name for outline colors
#'   (default \code{"viridis::viridis"}).
#' @param line_width Numeric line width for density outlines (default 1).
#' @return A ggplot object.
#' @export
plotZonationRidge <- function(zone_objs, palette = 'grDevices::rainbow', line_width = 1, adjust = 3) {
  if (is.list(zone_objs) && !is.null(zone_objs$mtx)) {
    zone_objs <- list(zone_objs)
  }
  if (is.null(names(zone_objs))) {
    names(zone_objs) <- paste0("Sample_", seq_along(zone_objs))
  } else {
    blank <- !nzchar(names(zone_objs))
    names(zone_objs)[blank] <- paste0("Sample_", which(blank))
  }
  df <- do.call(rbind, lapply(seq_along(zone_objs), function(i) {
    z <- as.numeric(getZonationGradient(zone_objs[[i]]))
    z <- z[is.finite(z)]
    if (length(z) == 0) {
      stop(sprintf("No finite zonation values for sample '%s'.",
                   names(zone_objs)[i]))
    }
    data.frame(Sample = names(zone_objs)[i], Zonation = z)
  }))
  df$Sample <- factor(df$Sample, levels = names(zone_objs))
  n <- length(zone_objs)
  pal_colors <- tryCatch({
    cols <- as.character(paletteer::paletteer_d(palette))
    rep_len(cols, n)
  }, error = function(e) {
    as.character(paletteer::paletteer_c(palette, n = n))
  })
  ggplot2::ggplot(df, ggplot2::aes(x = Zonation, color = Sample)) +
    ggplot2::geom_density(fill = NA, linewidth = line_width, bounds = c(1, 3), adjust = adjust) +
    ggplot2::scale_color_manual(values = pal_colors) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(1, 3)) +
    ggplot2::labs(x = "Zonation", y = "Proportion") +
    ggplot2::theme_gray()
}
