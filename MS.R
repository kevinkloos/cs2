MS <- function(test, model, type = "prob", manual_thres = 0.25) {
  # obtain rates
  if(type == "prob") {
    tpr_scores <- 
      model |>
      pluck("scores") |>
      filter(.true_class == 1) |>
      pull(prob_scores)
    
    fpr_scores <-
      model |>
      pluck("scores") |>
      filter(.true_class == 0) |>
      pull(prob_scores)
    
  }
  else if(type == "raw") {
    tpr_scores <- 
      model |>
      pluck("scores") |>
      filter(.true_class == 1) |>
      pull(raw_scores) 
    
    fpr_scores <-
      model |>
      pluck("scores") |>
      filter(.true_class == 0) |>
      pull(raw_scores)  
  }
  
  tpr <- function(x) {mean(tpr_scores >= x)}
  fpr <- function(x) {mean(fpr_scores >= x)}
  
  # helper function for adjusted count
  .ac <- function(cac, x) {
    out <- (cac - fpr(x)) / (tpr(x) - fpr(x))
  }
  
  # obtain test scores
  if(type == "raw") {
    .probabilities_to_distances <- function(mdl, probs) {
      prob_params <- mdl$fit$fit$fit@prob.model |> unlist()
      prob_to_dist <- 
        (log(1 - probs) - log(probs) - prob_params[2]) / prob_params[1]
      return(prob_to_dist)
    }
    
    if(model$fit$fit$fit$spec$method$fit$func['fun'] == "ksvm") {
      test_scores_pr <- map(test, ~ stats::predict(model$fit, .x, type = "prob") %>% pull(.pred_1), 
                         .progress = TRUE)
      test_scores <- map(test_scores_pr, ~ .probabilities_to_distances(model$fit, .x) %>% sort(),
                         .progress = TRUE)
    }
    else {
      test_scores <- map(test, ~ stats::predict(model$fit, .x, type = "raw") %>% sort(),
                         .progress = TRUE)
    }
  }
  else if(type == "prob") {
    test_scores <- 
      map(test, ~ stats::predict(model$fit, .x, type = "prob") %>% pull(.pred_1) %>% sort(),
          .progress = TRUE)
  }
  
  .get_median_sweep <- function(test_scores, manual_thres) {
    n_obs <- length(test_scores)
    all_cac <- seq(1, 1/n_obs, -1/n_obs) 
    filter_crit <- map_lgl(test_scores, ~ tpr(.x) - fpr(.x) >= manual_thres)
    
    incl_test <- subset(test_scores, filter_crit)
    incl_cac <- subset(all_cac, filter_crit)
    
    all_ac <- map2_dbl(incl_cac, incl_test, ~ .ac(.x, .y))
    out <- median(all_ac)
    return(out)
  }
  
  results <- map_dbl(test_scores, ~ .get_median_sweep(.x, manual_thres), .progress = TRUE)
  
  # truncate between 0 and 1
  out <-
    results %>%
    pmax(0) %>%
    pmin(1)
  
  return(out)
  
}

