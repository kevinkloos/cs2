library(tidyverse)


DyS <- function(test, model, n_bins = 8, tol = 10^(-5)) {
  # helper function to compute DyS
  .compute_histogram <- function(scores, n_bins, rng) {
    bns_selec <- cut(x = scores, breaks = seq(rng[1], rng[2], length.out = n_bins + 1))
    dns <- table(bns_selec) / length(bns_selec)
    names(dns) <- NULL
    dns <- as.vector(dns)
    return(dns)
  }
  
  .hellinger_distance <- function(p_vec, q_vec) {
    sqsq <- (sqrt(p_vec) - sqrt(q_vec))^2
    out <- sqrt(sum(sqsq)) / sqrt(2)
    return(out)
  }
  
  .tenary_search <- function(lft, rgt, h_dist, tol, h0, h1, h_test, mthd) {
    conv <- FALSE
    while(!conv) {
      lft3 <- lft + (rgt - lft) / 3
      rgt3 <- rgt - (rgt - lft) / 3
      
      if(h_dist(lft3, h0, h1, h_test, mthd) > h_dist(rgt3, h0, h1, h_test, mthd)) {
        lft <- lft3
      }
      else {
        rgt <- rgt3
      }
      conv <- abs(lft-rgt) < tol
    }
    out <- (lft + rgt) / 2
    return(out)
  }
  
  .histogram_distance <- function(prev_est, h0, h1, h_test, mthd) {
    h_train <- (1 - prev_est) * h0 + prev_est * h1
    dist_trte <- mthd(h_train, h_test)
    return(dist_trte)
  }
  
  .compute_dys <- function(tst_scr, h0, h1, n_bins, tol, rng, mthd) {
    h_test <- .compute_histogram(tst_scr, n_bins, rng)
    out <- .tenary_search(0, 1, .histogram_distance, tol, h0, h1, h_test, mthd)
    return(out)
  } 
  
  # obtain range and distance method, can be expanded to more
  rng <- c(0, 1)
  mthd <- get(".hellinger_distance")
  
  # obtain train and test scores (prob)
  pos_scores <- model |>
    pluck("scores") |>
    filter(.true_class == 1) |>
    pull(prob_scores) 
  
  neg_scores <- model |>
    pluck("scores") |>
    filter(.true_class == 0) |>
    pull(prob_scores) 
  
  test_scores <- 
    map(test, ~ stats::predict(model$fit, .x, type = "prob") |> pull(.pred_1) |> sort())
  
  # compute histogram of positive and negative scores
  h0 <- .compute_histogram(neg_scores, n_bins, rng)
  h1 <- .compute_histogram(pos_scores, n_bins, rng)
  
  # minimize distance mixture and test histograms with DyS
  out <- map_dbl(test_scores, ~ .compute_dys(.x, h0, h1, n_bins, tol, rng, mthd))
  return(out)
}
