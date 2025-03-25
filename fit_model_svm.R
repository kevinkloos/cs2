# fit_model_svm <- function(data_train, label_name, n_folds = 5, seed_nr = NULL,
#                           cost = tune::tune(), rbf_sigma = tune::tune()) {
#   
#   # Some stopping arguments
#   stopifnot("argument `data` should be a data frame." = "data.frame" %in% class(data_train))
#   stopifnot("argument `n_folds` should be less that the number of rows in your data." =
#               between(n_folds, 0, nrow(data_train)))
#   stopifnot("argument `n_folds` should be an integer." = n_folds %% 1 == 0)
#   stopifnot("argument `seed_nr` should be an integer." = seed_nr %% 1 == 0)
#   stopifnot("argument `seed_nr` should be at least one." = seed_nr >= 1)
#   
#   # Set seed number
#   if(!is.null(seed_nr)){set.seed(seed_nr)}
#   
#   # Split the labelled data into a training set and a test set
#   data_folds <- rsample::vfold_cv(data_train, v = n_folds, repeats = 1,
#                                   strata = !!enquo(label_name))
#   
#   # Define a basic linear SVM
#   svm_model <- 
#     parsnip::svm_rbf(mode = "classification", cost = cost, rbf_sigma = rbf_sigma) |> 
#     parsnip::set_engine("kernlab")
#   
#   # Obtain the formula
#   frm <- paste(label_name, "~ .") |> as.formula()
#   
#   # Add a recipe what we do with the data during the analysis
#   svm_recipe <- 
#     recipes::recipe(frm, data = data_train) |>
#     recipes::step_normalize()
#   
#   # Combine the model and the recipe in a workflow
#   svm_wflow <- 
#     workflows::workflow() %>% 
#     workflows::add_model(svm_model) |>
#     workflows::add_recipe(svm_recipe)
#   
#   # Fit the CV predictions and save them in a tibble
#   svm_resamples <-
#     svm_wflow |>
#     tune::fit_resamples(data_folds, control = tune::control_resamples(save_pred = TRUE))
#   
#   # Obtain the predicted scores of the cross-validated observations of each repeat
#   
#   pos_col <- 
#     paste0(".pred_", data_train |> pull(label_name) |> levels() |> unique() |> subset(c(F,T)))
#   
#   .probabilities_to_distances <- function(mdl, probs) {
#     prob_params <- mdl$fit$fit$fit@prob.model |> unlist()
#     prob_to_dist <- 
#       (log(1 - probs) - log(probs) - prob_params[2]) / prob_params[1]
#     return(prob_to_dist)
#   }
#   
#   final_model <- generics::fit(svm_wflow, data_train)
#   
#   # add the transformed scores (linear predictor, distance to decision boundary)
#   scores <-
#     svm_resamples |> 
#     dplyr::select(id, .predictions) |>
#     tidyr::unnest(cols = .predictions)  |>
#     dplyr::select(.true_class, prob_scores = as.name(pos_col)) |>
#     dplyr::mutate(raw_scores = .probabilities_to_distances(final_model, prob_scores))
#     
#   
#   return(list(fit = final_model,
#               scores = scores))
#   
# }

#############################


fit_model_svm <- function(data_train, label_name, n_folds = 5, seed_nr = NULL) {
  
  # Some stopping arguments
  stopifnot("argument `data` should be a data frame." = "data.frame" %in% class(data_train))
  stopifnot("argument `n_folds` should be less that the number of rows in your data." =
              between(n_folds, 0, nrow(data_train)))
  stopifnot("argument `n_folds` should be an integer." = n_folds %% 1 == 0)
  stopifnot("argument `seed_nr` should be an integer." = seed_nr %% 1 == 0)
  stopifnot("argument `seed_nr` should be at least one." = seed_nr >= 1)
  
  # Set seed number
  if(!is.null(seed_nr)){set.seed(seed_nr)}
  
  # Split the labelled data into a training set and a test set
  data_folds <- rsample::vfold_cv(data_train, v = n_folds, repeats = 1,
                                  strata = !!enquo(label_name))
  
  # Define a basic linear SVM
  svm_model <- 
    parsnip::svm_rbf(mode = "classification", cost = tune::tune(), rbf_sigma = tune::tune()) |> 
    parsnip::set_engine("kernlab")
  
  # Obtain the formula
  frm <- paste(label_name, "~ .") |> as.formula()
  
  # Add a recipe what we do with the data during the analysis
  svm_recipe <- 
    recipes::recipe(frm, data = data_train) |>
    recipes::step_normalize()
  
  # Combine the model and the recipe in a workflow
  svm_wflow <- 
    workflows::workflow() %>% 
    workflows::add_model(svm_model) |>
    workflows::add_recipe(svm_recipe)
  
  # Define a grid for hyperparameter tuning
  svm_grid <- 
    dials::grid_max_entropy(
      cost(), rbf_sigma(),
      size = 25
    )
  
  # Perform the grid search
  svm_tune <- tune::tune_grid(
    object = svm_wflow,
    resamples = data_folds,
    grid = svm_grid,
    metrics = yardstick::metric_set(roc_auc),
    control = tune::control_grid(verbose = TRUE)
  )
  
  # Extract the best parameters
  best_params <-
    svm_tune |> 
    tune::select_best(metric = "roc_auc")
  
  # And finalize the workflow with the best parameters
  svm_final <-
    svm_wflow |>
    tune::finalize_workflow(best_params)
  
  # Fit the CV predictions and save them in a tibble
  svm_resamples <-
    svm_final |>
    tune::fit_resamples(data_folds, control = tune::control_resamples(save_pred = TRUE))
  
  # Obtain the predicted scores of the cross-validated observations of each repeat
  pos_col <- 
    paste0(".pred_", data_train |> pull(label_name) |> levels() |> unique() |> subset(c(F,T)))
  
  .probabilities_to_distances <- function(mdl, probs) {
    prob_params <- mdl$fit$fit$fit@prob.model |> unlist()
    prob_to_dist <- 
      (log(1 - probs) - log(probs) - prob_params[2]) / prob_params[1]
    return(prob_to_dist)
  }
  
  # Fit the model with the best parameters
  final_model <- generics::fit(svm_final, data_train)
  
  # add the transformed scores (linear predictor, distance to decision boundary)
  scores <-
    svm_resamples |> 
    dplyr::select(id, .predictions) |>
    tidyr::unnest(cols = .predictions)  |>
    dplyr::select(.true_class, prob_scores = as.name(pos_col)) |>
    dplyr::mutate(raw_scores = .probabilities_to_distances(final_model, prob_scores))
  
  
  return(list(fit = final_model,
              scores = scores))
}
