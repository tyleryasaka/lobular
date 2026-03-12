em_zonation = function(mtx, init_w, iterations, density_cut, default_reg = 0.9, regularization = NULL, cor_thresh = 0, max_cells = 10000, verbose = FALSE) {
  if (!is.null(regularization)) {
    min_cor = mix_rate = rigidity = regularization
  } else {
    min_cor = mix_rate = rigidity = default_reg
  }
  if (ncol(mtx) > max_cells) {
    if (verbose) cat(sprintf("Subsampling %d / %d cells for training\n", max_cells, ncol(mtx)))
    idx = sample(ncol(mtx), max_cells)
    mtx = mtx[, idx]
  }
  mtx_sc = t(mtx)
  mtx_sc = scale(mtx_sc)
  mtx_sc[is.na(mtx_sc)] = 0
  n_cells = nrow(mtx_sc)
  n_genes = ncol(mtx_sc)
  n_sv = min(20, n_genes - 1)
  svd_res = irlba::irlba(mtx_sc, nv = n_sv)
  eigenvalues = svd_res$d^2 / (n_cells - 1)
  gamma = n_genes / n_cells
  mp_upper = (1 + sqrt(gamma))^2
  n_signals = sum(eigenvalues > mp_upper)
  if (verbose) {
    cat(sprintf("Parallel analysis (%d genes): %d signals detected (MP upper: %.2f)\n", n_genes, n_signals, mp_upper))
    cat(sprintf("Top 5 eigenvalues: %s\n", paste(round(eigenvalues[1:min(5, length(eigenvalues))], 2), collapse = ", ")))
  }
  signal_genes_filtered = NULL
  if (n_signals > 0) {
    V = svd_res$v[, 1:n_signals, drop = FALSE]
    rownames(V) = colnames(mtx_sc)
    max_loading = apply(V^2, 1, max)
    null_loading = 1 / n_genes
    init_in_mtx = intersect(names(init_w)[init_w != 0], colnames(mtx_sc))
    if (length(init_in_mtx) > 0) {
      init_loading = max_loading[init_in_mtx]
      if (verbose) {
        cat(sprintf("Init gene max PC loading - min: %.4f, median: %.4f, max: %.4f (null: %.4f)\n",
                    min(init_loading), median(init_loading), max(init_loading), null_loading))
        cat(sprintf("Init genes passing signal filter: %d / %d\n", sum(init_loading > null_loading * 2), length(init_in_mtx)))
      }
      signal_genes_filtered = names(init_loading)[init_loading > null_loading * 2]
    }
    common_sv = intersect(rownames(V), names(init_w))
    if (length(common_sv) > 0) {
      iw = init_w[common_sv]
      pc_cors = sapply(1:n_signals, function(k) {
        cor(V[common_sv, k], iw, method = "pearson")
      })
      if (verbose) {
        cat("PC correlations with init weights:\n")
        for (k in 1:n_signals) {
          cat(sprintf("  PC%d (eigenvalue %.2f): rho = %.3f\n", k, eigenvalues[k], pc_cors[k]))
        }
      }
      if (n_signals == 1 && is.null(regularization)) {
        min_cor = mix_rate = rigidity = 0.5
        if (verbose) cat("Only 1 signal detected, relaxing regularization.\n")
      }
    }
  }
  mtx_pos = t(mtx)
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
  if (cor_thresh > 0) {
    fr_unfilt = as.vector(mtx_pos %*% w)
    pt_unfilt = rank(fr_unfilt) / (length(fr_unfilt) + 1)
    gene_cor = apply(mtx_pos, 2, function(x) cor(x, pt_unfilt, method = "spearman"))
    sig_genes = names(gene_cor)[abs(gene_cor) > cor_thresh]
    w[!(names(w) %in% sig_genes)] = 0
    w = w / sqrt(sum(w^2))
    if (verbose) {
      cat(sprintf("Correlation filter (|rho| > %.2f): %d / %d genes retained\n",
                  cor_thresh, sum(w != 0), length(w)))
    }
  }
  fr = as.vector(mtx_pos %*% w)
  d = density(fr, n = 512)
  thresh = max(d$y) * density_cut
  lims_f = range(d$x[d$y > thresh])
  lo = lims_f[1]
  hi = lims_f[2]
  f_ecdf = ecdf(fr)
  g1 = names(w)[w < 0]
  g3 = names(w)[w > 0]
  get_sub_model = function(gs, w_vec, m_mtx) {
    ws = w_vec
    ws[!(names(w_vec) %in% gs)] = 0
    r = as.vector(m_mtx %*% ws)
    e_sub = ecdf(r)
    vals = e_sub(r)
    list(e = e_sub, e_min = min(vals), e_max = max(vals), genes = gs)
  }
  m1 = get_sub_model(g1, w, mtx_pos)
  m3 = get_sub_model(g3, w, mtx_pos)
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
    norm_provided = TRUE
  ))
}
