# Shared validation and preprocessing helpers --------------------------------
assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
}

# Estimate transformations from training data only.
transform_fit <- function(data, vars = VARS, log_vars = intersect(LOG_VARS, vars)) {
  x <- as.data.frame(data[, vars, drop = FALSE])
  for (v in log_vars) {
    if (any(x[[v]] <= 0, na.rm = TRUE)) stop("Non-positive value in log variable: ", v)
    x[[v]] <- log10(x[[v]])
  }
  center <- vapply(x, mean, numeric(1), na.rm = TRUE)
  scale_value <- vapply(x, sd, numeric(1), na.rm = TRUE)
  if (any(!is.finite(scale_value) | scale_value <= 0)) stop("Invalid scaling vector")
  list(vars = vars, log_vars = log_vars, center = center, scale = scale_value)
}

transform_apply <- function(data, fit) {
  x <- as.data.frame(data[, fit$vars, drop = FALSE])
  for (v in fit$log_vars) x[[v]] <- log10(x[[v]])
  out <- sweep(as.matrix(x), 2, fit$center, "-")
  out <- sweep(out, 2, fit$scale, "/")
  colnames(out) <- fit$vars
  out
}

# Projection, SOM diagnostics and partition helpers ---------------------------
nearest_centroid <- function(x, centers) {
  d2 <- outer(rowSums(x^2), rowSums(centers^2), "+") - 2 * tcrossprod(x, centers)
  max.col(-d2, ties.method = "first")
}

nearest_units <- function(x, codes) {
  d2 <- outer(rowSums(x^2), rowSums(codes^2), "+") - 2 * tcrossprod(x, codes)
  list(
    unit = max.col(-d2, ties.method = "first"),
    distance = sqrt(pmax(0, apply(d2, 1, min)))
  )
}

topographic_error <- function(model, x) {
  codes <- model$codes[[1]]
  d2 <- outer(rowSums(x^2), rowSums(codes^2), "+") - 2 * tcrossprod(x, codes)
  nearest <- t(apply(d2, 1, order))[, 1:2, drop = FALSE]
  unit_distance <- kohonen::unit.distances(model$grid)
  adjacent <- min(unit_distance[unit_distance > 0])
  mean(unit_distance[cbind(nearest[, 1], nearest[, 2])] > adjacent * 1.01)
}

neighbour_preservation <- function(model, x, neighbours = 10L) {
  n <- nrow(x)
  k <- min(neighbours, n - 1L)
  input_order <- apply(as.matrix(dist(x)), 1, order)
  unit_distance <- kohonen::unit.distances(model$grid)
  mapped <- model$unit.classif
  map_distance <- unit_distance[mapped, mapped, drop = FALSE]
  map_order <- apply(map_distance, 1, order)
  mean(vapply(seq_len(n), function(i) {
    a <- setdiff(input_order[, i], i)[seq_len(k)]
    b <- setdiff(map_order[, i], i)[seq_len(k)]
    length(intersect(a, b)) / k
  }, numeric(1)))
}

partition_som <- function(model, k, x_target = NULL) {
  code_groups <- cutree(hclust(dist(model$codes[[1]]), method = "ward.D2"), k = k)
  if (is.null(x_target)) return(code_groups[model$unit.classif])
  code_groups[nearest_units(x_target, model$codes[[1]])$unit]
}

# Agreement, consensus and uncertainty metrics -------------------------------
entropy <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (!length(p)) return(NA_real_)
  -sum(p * log(p))
}

adjusted_mutual_information <- function(x, y) {
  x <- as.integer(factor(x)); y <- as.integer(factor(y))
  tab <- table(x, y); n <- sum(tab)
  ai <- rowSums(tab); bj <- colSums(tab)
  mi <- 0
  for (i in seq_along(ai)) for (j in seq_along(bj)) {
    nij <- tab[i, j]
    if (nij > 0) mi <- mi + (nij / n) * log((n * nij) / (ai[i] * bj[j]))
  }
  hx <- -sum((ai / n) * log(ai / n)); hy <- -sum((bj / n) * log(bj / n))
  emi <- 0
  for (a in ai) for (b in bj) {
    lo <- max(1L, a + b - n); hi <- min(a, b)
    if (lo <= hi) for (nij in lo:hi) {
      prob <- exp(lchoose(a, nij) + lchoose(n - a, b - nij) - lchoose(n, b))
      emi <- emi + (nij / n) * log((n * nij) / (a * b)) * prob
    }
  }
  denom <- (hx + hy) / 2 - emi
  if (abs(denom) < .Machine$double.eps) return(1)
  (mi - emi) / denom
}

coassociation <- function(assignments) {
  n <- nrow(assignments)
  out <- matrix(0, n, n)
  valid <- matrix(0, n, n)
  for (j in seq_len(ncol(assignments))) {
    labels <- assignments[, j]
    ok <- !is.na(labels)
    out[ok, ok] <- out[ok, ok] + outer(labels[ok], labels[ok], "==")
    valid[ok, ok] <- valid[ok, ok] + 1
  }
  out / pmax(valid, 1)
}

consensus_partition <- function(assignments, k) {
  coassoc <- coassociation(assignments)
  distance <- 1 - coassoc
  distance[!is.finite(distance)] <- 1
  diag(distance) <- 0
  labels <- cutree(hclust(as.dist(distance), method = "average"), k = k)
  list(labels = labels, coassociation = coassoc)
}

sample_consensus_metrics <- function(labels, coassoc) {
  groups <- sort(unique(labels))
  rows <- lapply(seq_along(labels), function(i) {
    affinities <- vapply(groups, function(g) mean(coassoc[i, labels == g]), numeric(1))
    probabilities <- affinities / sum(affinities)
    own <- affinities[groups == labels[i]]
    other <- max(affinities[groups != labels[i]], na.rm = TRUE)
    data.frame(
      row_index = i,
      support = own,
      margin = own - other,
      entropy = entropy(probabilities) / log(length(groups))
    )
  })
  bind_rows(rows)
}

clusterwise_jaccard <- function(reference, candidate) {
  ref_groups <- sort(unique(reference))
  values <- vapply(ref_groups, function(g) {
    a <- reference == g
    max(vapply(sort(unique(candidate)), function(h) {
      b <- candidate == h
      sum(a & b) / sum(a | b)
    }, numeric(1)))
  }, numeric(1))
  data.frame(cluster = ref_groups, jaccard = values)
}

# Auxiliary internal-validity indices; these are not treated as ecological truth.
validity_metrics <- function(x, labels) {
  labels <- as.integer(factor(labels))
  k <- length(unique(labels)); n <- nrow(x)
  sil <- mean(cluster::silhouette(labels, dist(x))[, "sil_width"])
  center_all <- colMeans(x)
  within <- 0; between <- 0
  centers <- matrix(NA_real_, k, ncol(x))
  scatters <- numeric(k)
  for (g in seq_len(k)) {
    z <- x[labels == g, , drop = FALSE]
    centers[g, ] <- colMeans(z)
    within <- within + sum(rowSums((z - matrix(centers[g, ], nrow(z), ncol(z), byrow = TRUE))^2))
    between <- between + nrow(z) * sum((centers[g, ] - center_all)^2)
    scatters[g] <- mean(sqrt(rowSums((z - matrix(centers[g, ], nrow(z), ncol(z), byrow = TRUE))^2)))
  }
  ch <- (between / (k - 1)) / (within / (n - k))
  center_d <- as.matrix(dist(centers))
  ratios <- vapply(seq_len(k), function(i) {
    max(vapply(setdiff(seq_len(k), i), function(j) (scatters[i] + scatters[j]) / center_d[i, j], numeric(1)))
  }, numeric(1))
  data.frame(silhouette = sil, calinski_harabasz = ch, davies_bouldin = mean(ratios))
}

# Resample campaign-specific station blocks to preserve the observation hierarchy.
stratified_block_sample <- function(data, fraction, seed) {
  set.seed(seed)
  sampled <- data %>%
    distinct(campaign, block_id) %>%
    group_by(campaign) %>%
    slice_sample(prop = fraction, replace = FALSE) %>%
    ungroup()
  which(data$block_id %in% sampled$block_id)
}

# Fit each comparison model on the training matrix and assign target samples.
fit_kmeans_predict <- function(x_train, x_target, k, seed) {
  set.seed(seed)
  fit <- kmeans(x_train, centers = k, nstart = 50, iter.max = 200)
  list(labels = nearest_centroid(x_target, fit$centers), failure = NA_character_)
}

fit_ward_predict <- function(x_train, x_target, k) {
  train_labels <- cutree(hclust(dist(x_train), method = "ward.D2"), k = k)
  centers <- do.call(rbind, lapply(sort(unique(train_labels)), function(g) colMeans(x_train[train_labels == g, , drop = FALSE])))
  list(labels = nearest_centroid(x_target, centers), failure = NA_character_)
}

fit_gmm_predict <- function(x_train, x_target, k) {
  tryCatch({
    fit <- mclust::Mclust(x_train, G = k, verbose = FALSE)
    if (is.null(fit) || is.null(fit$parameters)) stop("No converged GMM")
    pred <- predict(fit, newdata = x_target)
    list(labels = pred$classification, failure = NA_character_, model_name = fit$modelName)
  }, error = function(e) list(labels = rep(NA_integer_, nrow(x_target)), failure = conditionMessage(e), model_name = NA_character_))
}

# Apply the prespecified clusterwise Jaccard interpretation thresholds.
stability_class <- function(jaccard) {
  cut(jaccard, breaks = c(-Inf, 0.60, 0.75, Inf), labels = c("fragile", "uncertain", "stable"), right = FALSE)
}
