em_zonation = function(mtx, init_w, iterations, density_cut, min_cor, mix_rate, rigidity = 1, norm_mtx = NULL, verbose = FALSE) {
  if (!is.null(norm_mtx)) {
    if (!identical(dim(norm_mtx), dim(mtx)))
      stop("norm_mtx must have the same dimensions as mtx")
    if (!identical(rownames(norm_mtx), rownames(mtx)))
      stop("norm_mtx must have the same row names (genes) as mtx")
    if (!identical(colnames(norm_mtx), colnames(mtx)))
      stop("norm_mtx must have the same column names (cells) as mtx")
  }
  rm = Matrix::rowMeans(mtx)
  raw_var = Matrix::rowMeans(mtx^2) - rm^2
  n_hvg = min(2000, nrow(mtx))
  hvg_idx = order(raw_var, decreasing = TRUE)[1:n_hvg]
  hvg_names = rownames(mtx)[hvg_idx]
  mtx_hvg = normalizeMatrix(mtx[hvg_idx, ])
  mtx_hvg = t(mtx_hvg)
  mtx_hvg = scale(mtx_hvg)
  mtx_hvg[is.na(mtx_hvg)] = 0
  n_cells = nrow(mtx_hvg)
  svd_res = svd(mtx_hvg, nu = 0, nv = min(20, ncol(mtx_hvg)))
  eigenvalues = svd_res$d^2 / (n_cells - 1)
  n_perms = 10
  n_compare = min(20, length(eigenvalues))
  null_eigs = matrix(0, n_perms, n_compare)
  for (p in 1:n_perms) {
    if (verbose) {
      cat(sprintf("Permutation %d/%d\n", p, n_perms))
    }
    mtx_perm = apply(mtx_hvg, 2, sample)
    sv_perm = svd(mtx_perm, nu = 0, nv = 0)$d
    null_eigs[p, ] = (sv_perm^2 / (n_cells - 1))[1:n_compare]
  }
  null_threshold = apply(null_eigs, 2, max)
  n_signals = 0
  for (j in 1:n_compare) {
    if (eigenvalues[j] > null_threshold[j]) n_signals = j
    else break
  }
  if (verbose) {
    cat(sprintf("Parallel analysis (%d HVGs): %d signals detected\n", n_hvg, n_signals))
    cat(sprintf("Top 5 eigenvalues: %s\n", paste(round(eigenvalues[1:min(5, length(eigenvalues))], 2), collapse = ", ")))
    cat(sprintf("Top 5 null thresholds: %s\n", paste(round(null_threshold[1:min(5, n_compare)], 2), collapse = ", ")))
  }
  signal_genes_filtered = NULL
  if (n_signals > 0) {
    V = svd_res$v[, 1:n_signals, drop = FALSE]
    rownames(V) = hvg_names
    max_loading = apply(V^2, 1, max)
    null_loading = 1 / n_hvg
    init_in_hvg = intersect(names(init_w)[init_w != 0], hvg_names)
    if (length(init_in_hvg) > 0) {
      init_loading = max_loading[init_in_hvg]
      if (verbose) {
        cat(sprintf("Init gene max PC loading - min: %.4f, median: %.4f, max: %.4f (null: %.4f)\n",
                    min(init_loading), median(init_loading), max(init_loading), null_loading))
      }
      signal_genes_filtered = names(init_loading)[init_loading > null_loading * 2]
      init_not_hvg = setdiff(names(init_w)[init_w != 0], hvg_names)
      if (verbose) {
        cat(sprintf("Init genes in HVGs: %d, not in HVGs: %d\n", length(init_in_hvg), length(init_not_hvg)))
        cat(sprintf("Init genes passing signal filter: %d / %d\n", length(signal_genes_filtered), length(init_in_hvg)))
      }
    }
  }
  if (!is.null(norm_mtx)) {
    mtx_pos = t(norm_mtx)
  } else {
    mtx_pos = t(normalizeMatrix(mtx))
  }
  g_m = colMeans(mtx_pos)
  mtx_pos = scale(mtx_pos, center = g_m, scale = FALSE)
  v = apply(mtx_pos, 2, var)
  keep = v > 0 & !is.na(v)
  mtx_pos = mtx_pos[, keep]
  g_m = g_m[keep]
  w = rep(0, ncol(mtx_pos))
  names(w) = colnames(mtx_pos)
  common = intersect(colnames(mtx_pos), names(init_w))
  w[common] = init_w[common]
  if (sum(w^2) > 0) w = w / sqrt(sum(w^2))
  init_w_kept = w
  if (!is.null(signal_genes_filtered) && length(signal_genes_filtered) >= 2) {
    init_w_kept[!(names(init_w_kept) %in% signal_genes_filtered)] = 0
    if (sum(init_w_kept^2) > 0) init_w_kept = init_w_kept / sqrt(sum(init_w_kept^2))
    w = init_w_kept
  }
  init_nonzero = names(init_w_kept)[init_w_kept != 0]
  for (i in 1:iterations) {
    g1 = names(w)[w < 0]
    g3 = names(w)[w > 0]
    w1 = w
    w1[!(names(w) %in% g1)] = 0
    w3 = w
    w3[!(names(w) %in% g3)] = 0
    r1 = as.vector(mtx_pos %*% w1)
    r3 = as.vector(mtx_pos %*% w3)
    pa = rank(rank(r1) + rank(r3)) / (nrow(mtx_pos) + 1)
    pc = pa - mean(pa)
    w_new = as.vector(t(mtx_pos) %*% pc)
    w_new = w_new / sqrt(sum(w_new^2, na.rm = TRUE))
    names(w_new) = colnames(mtx_pos)
    if (sum(w_new * w, na.rm = TRUE) < 0) w_new = -w_new
    w = w_new
    if (cor(w, init_w_kept, use = "complete.obs", method = "spearman") < min_cor) {
      w_mixed = (1 - mix_rate) * w + mix_rate * init_w_kept
      non_init = !(names(w_mixed) %in% init_nonzero)
      w_mixed[non_init] = w_mixed[non_init] * (1 - rigidity)
      w = w_mixed / sqrt(sum(w_mixed^2, na.rm = TRUE))
    }
  }
  fr = as.vector(mtx_pos %*% w)
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
  g1 = names(w)[w < 0]
  g3 = names(w)[w > 0]
  m1 = get_sub_model(g1, w, mtx_pos, fr, lo, hi)
  m3 = get_sub_model(g3, w, mtx_pos, fr, lo, hi)
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
    m3 = m3,
    norm_provided = is.null(norm_mtx)
  ))
}
