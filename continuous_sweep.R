continuous_sweep <- function(test, model, dens = "beta", type = "prob", 
                             manual_thres = NULL) {
  dens_options <- c("norm", "skew", "beta", "kernel")
  type_options <- c("raw", "prob")
  stopifnot("dens should contain: norm, skew, beta, kernel" = 
              dens %in% dens_options)
  stopifnot("type should contain: prob, raw" = type %in% type_options)
  stopifnot("choose one dens" = length(dens) == 1)
  stopifnot("choose one type" = length(type) == 1)
  if(is.numeric(manual_thres)) {
    stopifnot("manual_thres is either NULL or between zero and one." = 
                manual_thres >= 0 | manual_thres <= 1)
  }
  else {
    stopifnot("manual_thres is either NULL or between zero and one." = 
                is.null(manual_thres))
  }
  
  if((dens == "norm" | dens == "skew") & type == "prob") {
    stop("invalid combination between density and type of scores")
  }
  if(dens == "beta" & type == "raw") {
    stop("invalid combination between density and type of scores")
  }
  
  if("data.frame" %in% class(test)) {test <- list(test)}
  # helper functions
  
  ## cumulative density for skew-normal distribution
  pskew <- function(theta, mu, sdev, skew) {
    dlt <- skew / sqrt(1 + skew^2)
    b <- sqrt(pi / (pi - 2*dlt^2)) * sdev
    a <- mu - b * dlt * sqrt(2 / pi)  
    
    out <- 1 - sn::psn(x = theta, xi = a, omega = b, alpha = skew) %>% as.numeric()
    return(out)
  }
  
  pkernel <- function(theta, eval, est) {
    out <- approxfun(theta, eval, 1 - est, rule = 2)
    return(out)
  }
  
  ## estimate the parameters given a fitted model and a distribution
  estimate_params <- function(fitted_model, dens = "norm", type = "prob") {
    dens_options <- c("norm", "skew", "beta", "kernel")
    type_options <- c("raw", "prob")
    stopifnot("dens should contain: norm, skew, beta, kernel" = 
                dens %in% dens_options)
    stopifnot("type should contain: prob, raw" = type %in% type_options)
    stopifnot("choose one dens" = length(dens) == 1)
    stopifnot("choose one type" = length(type) == 1)
    
    if(type == "prob") {
      data_pos <-
        fitted_model %>%
        purrr::pluck("scores") %>%
        dplyr::filter(.true_class == 1) %>%
        dplyr::pull(prob_scores)
      
      data_neg <-
        fitted_model %>%
        purrr::pluck("scores") %>%
        dplyr::filter(.true_class == 0) %>%
        dplyr::pull(prob_scores)
    }
    else if(type ==  "raw") {
      data_pos <-
        fitted_model %>%
        purrr::pluck("scores") %>%
        dplyr::filter(.true_class == 1) %>%
        dplyr::pull(raw_scores)
      
      data_neg <-
        fitted_model %>%
        purrr::pluck("scores") %>%
        dplyr::filter(.true_class == 0) %>%
        dplyr::pull(raw_scores)
    }
    
    # possible options dependent on dens
    if (dens == "norm") {
      out <- list(
        pos = c(
          mean = mean(data_pos),
          sd = sd(data_pos),
          lower.tail = FALSE
        ),
        neg = c(
          mean = mean(data_neg),
          sd = sd(data_neg),
          lower.tail = FALSE
        ),
        rng = c(
          min(data_pos, data_neg),
          max(data_pos, data_neg)
        )
      )
      return(out)
    }
    if (dens == "beta") {
      param_pos <- 
        data_pos |> 
        EnvStats::ebeta() |> 
        pluck("parameters")
      
      param_neg <- 
        data_neg |> 
        EnvStats::ebeta() |> 
        pluck("parameters")
      
      
      out <- list(pos = c(param_pos, lower.tail = FALSE),
                  neg = c(param_neg, lower.tail = FALSE),
                  rng = c(min(data_pos, data_neg),
                          max(data_pos, data_neg)))

      return(out)
    }
    
    if (dens ==  "skew") {
      out <- list(
        pos = c(
          mu = mean(data_pos),
          sdev = sd(data_pos),
          skew = sn::msn.mle(y = data_pos)$dp$alpha
        ),
        neg = c(
          mu = mean(data_neg),
          sdev = sd(data_neg),
          skew = sn::msn.mle(y = data_neg)$dp$alpha
        ),
        rng = c(
          min(data_pos, data_neg),
          max(data_pos, data_neg)
        )
      )
      return(out)
    }
    
    if (dens == "kernel") {
      cd_pos <- ks::kcde(data_pos)
      cd_neg <- ks::kcde(data_neg)
      
      out <- list(pos = with(cd_pos, stats::approxfun(eval.points, 1 - estimate, rule = 2)),
                  neg = with(cd_neg, stats::approxfun(eval.points, 1 - estimate, rule = 2)),
                  rng = c(
                    min(data_pos, data_neg),
                    max(data_pos, data_neg)
                  ))
      
      return(out)
    }
  }
  
  ## compute adjusted count given a threshold, fitted parameters and distribution
  adjusted_count_cs <- function(cac, theta, train_scores, dens) {
    name <- paste0("p", dens, collapse = "")
    if(name == "pkernel") {
      tpr <- train_scores$pos
      fpr <- train_scores$neg
    }
    else {
      tpr <- function(theta) {map_dbl(theta, ~ do.call(name, as.list(c(.x, train_scores$pos))))}
      fpr <- function(theta) {map_dbl(theta, ~ do.call(name, as.list(c(.x, train_scores$neg))))}
    }
    
    tpr_val <- map_dbl(theta, ~ tpr(.x))
    fpr_val <- map_dbl(theta, ~ fpr(.x))
    
    return((cac - fpr_val) / (tpr_val - fpr_val))
  }
  
  ## calculate the area under the curve between two thresholds
  calculate_one_area <- function(cac, theta_l, theta_r, train_scores, dens) {
    
    out <- stats::integrate(f = function(x) {
      adjusted_count_cs(cac, x, train_scores, dens)
    }, lower = theta_l, upper = theta_r)$value
    
    return(out)
  }
  
  ## get the thresholds given the fitted parameters, distribution, and cut-off value
  get_threshold_barriers <- function(train_scores, dens, pdelta) {
    name <- paste0("p", dens, collapse = "")
    if(name == "pkernel") {
      tpr <- train_scores$pos
      fpr <- train_scores$neg
    }
    else {
      tpr <- function(theta) {map_dbl(theta, ~ do.call(name, as.list(c(.x, train_scores$pos))))}
      fpr <- function(theta) {map_dbl(theta, ~ do.call(name, as.list(c(.x, train_scores$neg))))}
    }
    diff_rates <- function(theta) {tpr(theta) - fpr(theta)}
    
    barriers <-
      rootSolve::uniroot.all(f = function(x) {
        out <- tpr(x) - fpr(x) - pdelta
        return(out)
      }, interval = train_scores$rng)
    
    if(length(barriers) == 0) {
      lft_range <- tpr(train_scores$rng[1]) - fpr(train_scores$rng[1])
      rgt_range <- tpr(train_scores$rng[2]) - fpr(train_scores$rng[2])
      if(lft_range >= pdelta && rgt_range >= pdelta) {
        barriers <- train_scores$rng
      }
      else {
        return(NA)
      }
    }
    
    if(length(barriers) %% 2 == 1) {
      lft <- train_scores$rng[1]
      rgt <- train_scores$rng[2]
      # append it with otherwise the left or right boundary 
      # (whichever is suitable)
      if(tpr(lft) - fpr(lft) > pdelta) {
        barriers <- c(lft, barriers)
      }
      else if(tpr(rgt) - fpr(rgt) > pdelta){
        barriers <- c(barriers, rgt)
      }
    }
    return(barriers)
  }
  
  ## integrate to compute the covariance of continuous sweep
  continuous_sweep_int <- function(n_obs, train_scores, dens, trunc_value) {
    # obtain decision boundaries
    barrier_levels <- get_threshold_barriers(train_scores, dens, trunc_value)
    # compute integral
    intgrl <- 
      pracma::integral2(fun = function(xi, yi) {
        adjusted_count_cov(xi, yi, n_obs, train_scores, dens)
      }, barrier_levels[1], barrier_levels[2], barrier_levels[1], barrier_levels[2])$Q 
    
    scaling_factor <- diff(barrier_levels)^2
    out <- intgrl / scaling_factor
    
    return(out)
  }
  
  ## obtain maximum difference between tpr and fpr (needed for optimal threshold)
  get_optim_threshold <- function(train_scores, dens) {
    name <- paste0("p", dens, collapse = "")
    if(name == "pkernel") {
      tpr <- train_scores$pos
      fpr <- train_scores$neg
    }
    else {
      tpr <- function(theta) {map_dbl(theta, ~ do.call(name, as.list(c(.x, train_scores$pos))))}
      fpr <- function(theta) {map_dbl(theta, ~ do.call(name, as.list(c(.x, train_scores$neg))))}
    }
    rng <- train_scores$rng
    diff_rates <- function(theta) {tpr(theta) - fpr(theta)}
    
    # avoid boundaries
    err <- 10^(-5)
    
    max_pdelta <- 
      optimize(f = function(x) {tpr(x) - fpr(x)}, 
               interval = rng, maximum = TRUE) |>
      purrr::pluck("objective")
    
    min_pdelta <- 
      pmax(tpr(rng[1]) - fpr(rng[1]),
           tpr(rng[2]) - fpr(rng[2]))
    
    .estimate_mse <- function(tpr, fpr, rng, pdelta) {
      # compute covariance at point x, y
      .cov_ac <- function(x_coor, y_coor, rates_cs) {
        prev <- 0.5
        diff_rates <- function(x) {tpr(x) - fpr(x)}
        
        pt1 <- purrr::map2_dbl(x_coor, y_coor, ~ tpr(max(.x, .y)) - tpr(.x) * tpr(.y))
        pt2 <- purrr::map2_dbl(x_coor, y_coor, ~ fpr(max(.x, .y)) - fpr(.x) * fpr(.y))
        
        cov_th <-  prev * (pt1 + pt2)
        
        diff_rates_x <- purrr::map_dbl(x_coor, ~ diff_rates(.x))
        diff_rates_y <- purrr::map_dbl(y_coor, ~ diff_rates(.x))
        
        out <- cov_th / (diff_rates_x * diff_rates_y) 
        return(out)
      }
      # find the roots where F^+ - F^- = pdelta
      roots_pdelta <- 
        rootSolve::uniroot.all(f = function(x) {
          out <- tpr(x) - fpr(x) - pdelta
          return(out)
        }, interval = rng)
      
      # check whether length roots pdelta is not zero, otherwise infinite variance
      if(length(roots_pdelta) == 0) {
        return(Inf)
      }
      
      # check whether roots pdelta is an odd number
      if(length(roots_pdelta) %% 2 == 1) {
        # append it with otherwise the left or right boundary (whichever is suitable)
        if(tpr(rng[1]) - fpr(rng[1]) > pdelta) {
          roots_pdelta <- c(rng[1], roots_pdelta)
        }
        else if(tpr(rng[2]) - fpr(rng[2]) > pdelta){
          roots_pdelta <- c(roots_pdelta, rng[2])
        }
      }
      # sum the integrals (also works if there is one interval)
      cov_auc <- numeric(length(roots_pdelta) / 2)
      idx_sum <- 1
      for(idx in seq(1, length(roots_pdelta), 2)) {
        theta_l <- roots_pdelta[idx]
        theta_r <- roots_pdelta[idx + 1]
        if(theta_l == theta_r) {
          cov_auc[idx_sum] <- .cov_ac(theta_l, theta_r, rates_cs)
        }
        else{
          cov_auc[idx_sum] <- pracma::integral2(f = function(x, y) {
            .cov_ac(x, y, rates_cs)
          }, theta_l, theta_r, theta_l, theta_r)$Q / (theta_r - theta_l)^2
        }
        idx_sum <- idx_sum + 1
      }
      
      return(sum(cov_auc))
    }
    
    # need estimate of mse
    optim_pdelta <- 
      optimise(f = function(x) {.estimate_mse(tpr, fpr, rng, x)},
               interval = c(min_pdelta + err, max_pdelta - err)) |>
      purrr::pluck("minimum")
    
    return(optim_pdelta)
  }
  
  ## compute continuous sweep
  
  get_continuous_sweep <- function(train_scores, test_scores, dens, trunc_value) {
    # set variables
    n_obs <- length(test_scores)
    
    # make 3 vectors containing the left bound, right bound, and cac-value
    all_theta_l <- test_scores %>% utils::head(-1)
    all_theta_r <- test_scores %>% utils::tail(-1)
    all_cac <- seq(1 - 1/n_obs, 1/n_obs, -1/n_obs) 
    
    barriers <- get_threshold_barriers(train_scores, dens, trunc_value)
    area_sum <- numeric(length(barriers) / 2)
    diff_range <- numeric(length(barriers) / 2)
    
    # filter out the intervals that fall outside the boundaries
    # suitable for more than a pair of decision boundaries
    for(idx in seq(from = 1, to = length(barriers), by = 2)) {
      idx_sum <- 1
      diff_range[idx_sum] <- barriers[idx + 1] - barriers[idx]
      bool_crit <- 
        all_theta_r >= barriers[idx] & all_theta_l <= barriers[idx + 1]
      
      if(all(bool_crit == FALSE)){
        area_sum[idx_sum] <- NULL
        diff_range[idx_sum] <- NULL
        idx_sum <- idx_sum + 1
        next
      }
      n_trunc <- sum(bool_crit)
      
      # replace the values slightly below and sightly above barriers 
      # with correct values
      all_theta_l <-
        subset(all_theta_l, bool_crit) |> 
        replace(1, barriers[idx])
      all_theta_r <-
        subset(all_theta_r, bool_crit) |> 
        replace(n_trunc, barriers[idx + 1])
      all_cac <- subset(all_cac, bool_crit)
      
      all_areas <- 
        purrr::pmap_dbl(.l = list(all_cac, all_theta_l, all_theta_r), 
                        .f =  ~ calculate_one_area(..1, ..2, ..3, 
                                                   train_scores, dens))
      
      area_sum[idx_sum] <- sum(all_areas)
      idx_sum <- idx_sum + 1
    }
    if(all(is.null(area_sum))){return(NA)}
    return(sum(area_sum) / sum(diff_range))
  }
  
  # scores
  train_scores <- estimate_params(model, dens, type)
  
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
      test_scores <- map(test, ~ stats::predict(model$fit, .x, type = "raw") %>% sort())
    }
  }
  else if(type == "prob") {
    test_scores <- 
      map(test, ~ stats::predict(model$fit, .x, type = "prob") %>% pull(.pred_1) %>% sort())
  }
  
  # obtain thresholds
  optim_thres <- get_optim_threshold(train_scores, dens)
  
  # estimate prevalence for each test set
  if(!is.null(manual_thres)) {
    results <-
      map_dbl(test_scores,
              ~ get_continuous_sweep(train_scores, .x, dens, manual_thres))
  }
  else {
    results <-
      map_dbl(test_scores,
              ~ get_continuous_sweep(train_scores, .x, dens, optim_thres))
  }
  
  # truncate between 0 and 1
  out <-
    results %>%
    pmax(0) %>%
    pmin(1)
  
  return(out)
}
