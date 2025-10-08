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
ZonationObject = function(baseline, species, factors) {
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

  # Create the object
  obj = list(
    baseline = baseline,
    species = species,
    factors = factors
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
  invisible(object)
}

#' Convert mouse genes to human
#'
#' @param mouse_genes A vector of mouse gene names
#' @return A vector of human gene names
#' @noRd
mouseToHuman = function(mouse_genes) {
  mouse_entrez = AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                        keys = mouse_genes,
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first")
  human_genes_simple = toupper(mouse_genes)
  human_genes_verified = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                keys = human_genes_simple,
                                column = "SYMBOL",
                                keytype = "SYMBOL",
                                multiVals = "first")
  result = data.frame(
    mouse_symbol = mouse_genes,
    mouse_entrez = mouse_entrez,
    human_symbol = human_genes_verified,
    stringsAsFactors = FALSE
  )
  return(result$human_symbol)
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

#' Calibrate the model to baseline liver zonation
#' @param mtx Gene expression matrix with genes as rows
#' @param species Species to use, defaults to human. Currently supports 'mouse' and 'human'.
#' @return A \code{ZonationObject} with calibrated baseline zonation
#' @export
setBaseline = function(mtx, species = 'human') {
  hep_zonated = read.csv(system.file('extdata', 'xu_naturegenetics_2024_hepatocyte_zonated_genes.csv', package = 'lobular'))
  layers = sapply(1:9, function(x) paste0('Layer.', x))
  hep_zonated$zonation = apply(hep_zonated, 1, function(row) {
    vec = as.numeric(row[layers])
    vec = expanded_data = rep(1:length(vec), times = round(vec * 100))
    moments::skewness(vec)
  })

  if (species == 'human') {
    hep_zonated$Gene_converted = mouseToHuman(hep_zonated$Gene)
  } else if (species == 'mouse') {
    hep_zonated$Gene_converted = hep_zonated$Gene
  } else {
    stop("Only 'human' and 'mouse' species are supported at the moment. (Specify with species = 'mouse'")
  }
  genes.zonated = intersect(hep_zonated$Gene_converted, rownames(mtx))
  hep_zonated = hep_zonated[hep_zonated$Gene_converted %in% genes.zonated,]
  factors.zonated = hep_zonated$zonation
  names(factors.zonated) = hep_zonated$Gene_converted

  zonescore_raw = getGeneAvg(mtx, factors.zonated)
  zone_obj = ZonationObject(
    baseline = zonescore_raw,
    species = species,
    factors = factors.zonated
  )
  zone_obj
}

#' Apply the model to new values, returning the zonation
#' @param mtx Gene expression matrix with genes as rows
#' @param zone_obj Calibrated Zonation Object
#' @return A vector of zonation assignments
#' @export
getZone = function(mtx, zone_obj) {
  new_vec = getGeneAvg(mtx, zone_obj$factors)
  zonescore = apply_transformation(new_vec, zone_obj$baseline)
  zone = ifelse(zonescore < (1 / 3), 1, ifelse(zonescore < (2 / 3), 2, 3))
  zone = factor(zone, levels = c('1', '2', '3'))
  zone
}
