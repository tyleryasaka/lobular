ZONES = c('Zone_1', 'Zone_2', 'Zone_3')
ZONE_COLORS = c('#504880FF', '#1BB6AFFF', '#FFAD0AFF', '#D72000FF', '#F080D8FF')
minMaxNorm = function(v) (v - min(v, na.rm = T)) / (max(v, na.rm = T) - min(v, na.rm = T))

em_zonation = function(mtx, init_w, iterations, density_cut, min_cor, mix_rate) {
  mtx = t(mtx)
  g_m = colMeans(mtx)
  mtx = scale(mtx, center = g_m, scale = FALSE)
  v = apply(mtx, 2, var)
  keep = v > 0 & !is.na(v)
  mtx = mtx[, keep]
  g_m = g_m[keep]
  w = rep(0, ncol(mtx))
  names(w) = colnames(mtx)
  common = intersect(colnames(mtx), names(init_w))
  w[common] = init_w[common]
  if (sum(w^2) > 0) w = w / sqrt(sum(w^2))
  init_w_kept = w
  for (i in 1:iterations) {
    g1 = names(w)[w < 0]
    g3 = names(w)[w > 0]
    w1 = w
    w1[!(names(w) %in% g1)] = 0
    w3 = w
    w3[!(names(w) %in% g3)] = 0
    r1 = as.vector(mtx %*% w1)
    r3 = as.vector(mtx %*% w3)
    pa = rank(rank(r1) + rank(r3)) / (nrow(mtx) + 1)
    pc = pa - mean(pa)
    w_new = as.vector(t(mtx) %*% pc)
    w_new = w_new / sqrt(sum(w_new^2, na.rm = TRUE))
    names(w_new) = colnames(mtx)
    if (sum(w_new * w, na.rm = TRUE) < 0) w_new = -w_new
    w = w_new
    if (cor(w, init_w_kept, use = "complete.obs", method = "spearman") < min_cor) {
      w = (1 - mix_rate) * w + mix_rate * init_w_kept
      w = w / sqrt(sum(w^2, na.rm = TRUE))
    }
  }
  fr = as.vector(mtx %*% w)
  d = density(fr, n = 512)
  thresh = max(d$y) * density_cut
  lims_f = range(d$x[d$y > thresh])
  lo = lims_f[1]
  hi = lims_f[2]
  f_ecdf = ecdf(fr)
  get_sub_model = function(gs, w_vec, m_mtx, f_vec, l_val, h_val) {
    ws = w_vec
    ws[!(names(w_vec) %in% gs)] = 0
    r = as.vector(m_mtx %*% ws)
    b = coef(lm(r ~ f_vec))
    e_sub = ecdf(r)
    list(e = e_sub, rl = e_sub(b[1] + b[2] * l_val), rh = e_sub(b[1] + b[2] * h_val))
  }
  m1 = get_sub_model(g1, w, mtx, fr, lo, hi)
  m3 = get_sub_model(g3, w, mtx, fr, lo, hi)
  pt = (f_ecdf(fr) - f_ecdf(lo)) / (f_ecdf(hi) - f_ecdf(lo))
  return(list(
    weights = w,
    gene_means = g_m,
    f_e = f_ecdf,
    f_rl = f_ecdf(lo),
    f_rh = f_ecdf(hi),
    q1 = quantile(pt, 1/3, na.rm = TRUE),
    q2 = quantile(pt, 2/3, na.rm = TRUE),
    m1 = m1,
    m3 = m3
  ))
}

predict_position = function(mtx, model) {
  common_genes = intersect(rownames(mtx), names(model$weights))
  mtx = mtx[common_genes, , drop = FALSE]
  mtx = scale(t(mtx), center = model$gene_means[common_genes], scale = FALSE)
  w_s = model$weights[common_genes]
  if(length(w_s) > 0) w_s = w_s / sqrt(sum(w_s^2))
  raw = as.vector(mtx %*% w_s)
  (model$f_e(raw) - model$f_rl) / (model$f_rh - model$f_rl)
}

predict_position_2d = function(mtx, model) {
  g1 = names(model$weights)[model$weights < 0]
  g3 = names(model$weights)[model$weights > 0]
  pg = predict_position(mtx, model)
  common_genes = intersect(rownames(mtx), names(model$weights))
  mtx_c = scale(t(mtx[common_genes, , drop = FALSE]), center = model$gene_means[common_genes], scale = FALSE)
  gp = function(gs, sub_m) {
    valid_gs = intersect(gs, common_genes)
    w_s = model$weights[common_genes]; w_s[!(names(w_s) %in% valid_gs)] = 0
    raw = as.vector(mtx_c %*% w_s)
    (sub_m$e(raw) - sub_m$rl) / (sub_m$rh - sub_m$rl)
  }
  z = ifelse(pg > model$q2, 'Zone_3', ifelse(pg > model$q1, 'Zone_2', 'Zone_1'))
  data.frame(ZONE_1 = 1 - gp(g1, model$m1), ZONE_3 = gp(g3, model$m3),
             zonation = pg, zone = factor(z, levels = c('Zone_1', 'Zone_2', 'Zone_3')))
}




normalizeMatrix = function(mtx) {
  as.matrix(mtx) # no normalization (user should pass in log-transformed matrix)
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

getZonationGradient_help = function(mtx, zone_obj) {
  predict_position(mtx, zone_obj) * 2 + 1
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
setBaseline = function(mtx, coords = NULL, species = 'human', min_cor = 0.5, mix_rate = 0.5) {
  if (species == 'human') {
    initial_weights = readRDS(system.file('extdata', 'initial_weights_human.RDS', package = 'lobular'))
  } else if (species == 'mouse') {
    initial_weights = readRDS(system.file('extdata', 'initial_weights_mouse.RDS', package = 'lobular'))
  } else {
    stop("Only 'human' and 'mouse' species are supported at the moment. (Specify with species = 'mouse'")
  }
  # initial_weights = initial_weights[abs(initial_weights) > 0.025]
  mtx = as.matrix(mtx)
  em_zonation(mtx, initial_weights, iterations = 10, density_cut = 0, min_cor = min_cor, mix_rate = mix_rate)
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
  mtx = normalizeMatrix(mtx)
  getZonationGradient_help(mtx, zone_obj)
}

#' Apply the model to new samples, returning the zonation per cell/spot as discrete bins (1, 2, or 3)
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A vector of zonation assignments (discrete)
#' @export
getZone = function(mtx, zone_obj) {
  mtx = normalizeMatrix(mtx)
  zonescore = predict_position(mtx, zone_obj)
  zone = ifelse(zonescore < (1 / 3), 'Zone_1', ifelse(zonescore < (2 / 3), 'Zone_2', 'Zone_3'))
  zone = factor(zone, levels = ZONES)
  names(zone) = colnames(mtx)
  zone
}

getZonation2d_help = function(mtx, zone_obj) {
  zone.2d = predict_position_2d(mtx, zone_obj)
  rownames(zone.2d) = colnames(mtx)
  zone.2d
}

#' Apply the model to new samples, returning the 2d zonation per cell/spot
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A dataframe with columns PV, CV, and zone
#' @export
getZonation2d = function(mtx, zone_obj) {
  mtx = normalizeMatrix(mtx)
  getZonation2d_help(mtx, zone_obj)
}

#' Apply the model to new samples, returning a plot of the 2d zonation per cell/spot with *zone* indicated by the color
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @param point_size Optional numeric value for the ggplot point size (default 1)
#' @return A ggplot object
#' @export
plotZonation2d = function(mtx, zone_obj, point_size = 1) {
  mtx = normalizeMatrix(mtx)
  zone.2d = getZonation2d_help(mtx, zone_obj)
  ggplot(zone.2d) + geom_point(aes(ZONE_1, ZONE_3, color = zone), size = point_size) + coord_fixed() +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
    geom_abline(slope = 1, intercept = 1/3, color = '#D72000FF') +
    geom_abline(slope = 1, intercept = -1/3, color = '#1BB6AFFF') +
    geom_abline(slope = -1, intercept = 1/3, color = '#504880FF') +
    geom_abline(slope = -1, intercept = 5/3, color = '#F080D8FF') +
    scale_color_manual(values = ZONE_COLORS, limits = ZONES) +
    labs(x = 'Zone 1 Score', y = 'Zone 3 Score', color = 'Zone') +
    ggdark::dark_theme_grey()
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
  mtx = normalizeMatrix(mtx)
  zone.2d = getZonation2d_help(mtx, zone_obj)
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
plotZonation2dGene = function(mtx, zone_obj, gene, point_size = 1) {
  mtx = normalizeMatrix(mtx)
  zone.2d = getZonation2d_help(mtx, zone_obj)
  zone.2d[,gene] = mtx[gene,]
  ggplot(zone.2d) + geom_point(aes(ZONE_1, ZONE_3, color = .data[[gene]]), size = point_size) + coord_fixed() +
    scale_x_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
    scale_color_viridis_c() +
    labs(x = 'Zone 1 Score', y = 'Zone 3 Score', color = gene) +
    ggdark::dark_theme_grey()
}

#' Apply the model to new samples, returning a ridge plot of gene expression along zonation axis
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @param gene Name of gene to plot
#' @return A ggplot object
#' @export
plotRegression = function(mtx, zone_obj, gene) {
  mtx = normalizeMatrix(mtx)
  zone.2d = getZonation2d_help(mtx, zone_obj)
  zone.2d[,gene] = mtx[gene,]
  ggplot(zone.2d, aes(x = zonation, y = .data[[gene]])) +
    geom_point(color = '#504880FF') +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = TRUE, color = "#F080D8FF", fill  = "#C8C0F8FF") +
    xlab('Zonation') +
    ggdark::dark_theme_classic()
}

#' Apply the model to new samples, returning a ridge plot of the scores for zones 1, 2, and 3
#'
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A ggplot object
#' @export
plotZonationRidge = function(mtx, zone_obj) {
  mtx = normalizeMatrix(mtx)
  zone.2d = getZonation2d_help(mtx, zone_obj)
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
  mtx = normalizeMatrix(mtx)
  zone.2d = getZonation2d_help(mtx, zone_obj)
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
getZoneSpatial = function(mtx, coords, zone_obj, resolution = 1, use_for_inference = NULL) {
  mtx = normalizeMatrix(mtx)
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
plotZoneSpatial = function(mtx, coords, zone_obj, resolution = 1, use_for_inference = NULL) {
  mtx = normalizeMatrix(mtx)
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
  mtx = normalizeMatrix(mtx)
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
  mtx = normalizeMatrix(mtx)
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
