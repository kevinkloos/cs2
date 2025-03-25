SLD <- function(test, model) {
  # compute scores
  alpha_train <- 
    model %>%
    pluck("scores") %>%
    count(.true_class) %>%
    mutate(probs = n / sum(n)) %>%
    pull(probs)
  
  test_scores <- 
    map(test, ~ stats::predict(model$fit, .x, type = "prob") %>% pull(.pred_1) %>% sort())
  
  SLD_alg <- function(prob_test, alpha_train, maxiter = 1000, crit = 0.0001) { 
    # maxiter and crit from QuaPy

    prob_test <- cbind(1 - prob_test, prob_test)
    alpha_train <- alpha_train %>% matrix(nrow = 1)
    alpha_new <- alpha_train
    prob_test <- as.matrix(prob_test)
    iter <- 0
    diff <- Inf
    
    .updater <- function(alpha_new, alpha_train, prob_test) {
      n <- length(alpha_train)
      # E-step
      unscaled <- 
        map(1:2, ~ alpha_new[.x] / alpha_train[.x] * prob_test[,.x]) %>%
        bind_cols(.name_repair = "unique") %>%
        suppressMessages()
      scaled <- unscaled / rowSums(unscaled)
      
      # U-step
      out <- colMeans(scaled)
      names(out) <- NULL
      return(out)
    }
    
    while(iter < maxiter && crit <= diff) {
      alpha_upd <- .updater(alpha_new, alpha_train, prob_test)
      diff <- abs(alpha_upd - alpha_new) %>% mean()
      iter <- iter + 1
      alpha_new <- alpha_upd
    }
    return(alpha_new)
  }
  
  test_predict_SLD <- 
    map_dbl(test_scores, ~ SLD_alg(.x, alpha_train) %>% dplyr::last(), 
            .progress = TRUE)
  
  return(test_predict_SLD)
}
