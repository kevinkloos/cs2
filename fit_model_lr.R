fit_model_lr <- function(data, label_name, n_folds = 5,
                         seed_nr = NULL, penalty = 0, mixture = 1) {
  # Some stopping arguments
  stopifnot("argument `data` should be a data frame." = "data.frame" %in% class(data))
  stopifnot("argument `n_folds` should be less that the number of rows in your data." =
            between(n_folds, 0, nrow(data)))
  stopifnot("argument `n_folds` should be an integer." = n_folds %% 1 == 0)
  stopifnot("argument `seed_nr` should be an integer." = seed_nr %% 1 == 0)
  stopifnot("argument `seed_nr` should be at least one." = seed_nr >= 1)
  stopifnot("argument `penalty` should be at least zero." = penalty >= 0)
  stopifnot("argument `mixture` shuld be between zero and one." = between(mixture, 0, 1))
  
  # Set seed number
  if(!is.null(seed_nr)){set.seed(seed_nr)}
  
  # Split the labelled data into a training set and a test set
  data_folds <- rsample::vfold_cv(data, v = n_folds, strata = {{ label_name }})
  
  # Define a basic logistic regression model
  # check penalty (Leuven!!!)
  lr_model <- 
    parsnip::logistic_reg(penalty = penalty, mixture = mixture) %>% 
    parsnip::set_engine("glmnet")
  
  # Obtain the formula
  frm <- paste(deparse(substitute(label_name)), "~ .") %>% as.formula()
  
  # Add a recipe what we do with the data during the analysis
  lr_recipe <- 
    recipes::recipe(frm, data = data)
  
  # Combine the model and the recipe in a workflow
  lr_wflow <- 
    workflows::workflow() %>% 
    workflows::add_model(lr_model) %>% 
    workflows::add_recipe(lr_recipe)
  
  # Fit the best model
  lr_fit <- 
    lr_wflow %>%
    generics::fit(data = data)
  
  lr_resamples <- 
    lr_wflow %>%
    tune::fit_resamples(data_folds, control = tune::control_resamples(save_pred = TRUE)) %>%
    tune::collect_predictions()
  
  # Return the scores of the positive and negative class
  scores <- 
    lr_resamples %>%
    dplyr::select(.pred_1, .pred_class, .true_class) %>%
    dplyr::rename(prob_scores = .pred_1) %>%
    dplyr::mutate(raw_scores = log(prob_scores / (1 - prob_scores)), .after = prob_scores)
    
  return(list(scores = scores, fit = lr_fit))
}
