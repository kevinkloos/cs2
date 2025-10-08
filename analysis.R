source("fit_model_lr.R")
source("fit_model_svm.R")
source("continuous_sweep.R")
source("SLD.R")
source("MS.R")
source("DyS.R")
library(tidyverse)
library(tidymodels)

train <- 
  readr::read_delim("dat/training_data.txt", show_col_types = FALSE) %>%
  dplyr::rename_at(.vars = vars(matches('[0-9]')), .funs = ~ paste0('V', .)) %>%
  dplyr::rename(.true_class = label) %>%
  dplyr::mutate(.true_class = as.factor(.true_class))

test <- 
  map(0:4999, function(.x) {
    readr::read_delim(paste0("dat/test_samples/", .x, ".txt"), show_col_types = FALSE) %>%
      dplyr::rename_at(.vars = vars(matches('[0-9]')), .funs = ~ paste0('V', .))
  }, .progress = TRUE)

# set.seed(1234)

# fitted_lr <- fit_model_lr(train, .true_class, penalty = 0.006, seed_nr = 123)
# fitted_svm <- fit_model_svm(train, ".true_class", seed_nr = 123)
fitted_svm <- readRDS("fitted_svm.RData")

#opt thres
prev_norm <- continuous_sweep(test, fitted_svm, "norm", "raw")
cat("prev_norm completed \n")
timestamp()
prev_skew <- continuous_sweep(test, fitted_svm, "skew", "raw")
cat("prev_skew completed \n")
timestamp()
prev_kernelr <- continuous_sweep(test, fitted_svm, "kernel", "raw")
cat("prev_kernelr completed \n")
timestamp()
prev_kernelp <- continuous_sweep(test, fitted_svm, "kernel", "prob")
cat("prev_kernelp completed \n")
timestamp()
prev_beta <- continuous_sweep(test, fitted_svm, "beta", "prob")
cat("prev_beta completed \n")
timestamp()

#trad thres
prev_norm_trad <- continuous_sweep(test, fitted_svm, "norm", "raw",
                                   manual_thres = 0.25)
cat("prev_norm_trad completed \n")
timestamp()
prev_skew_trad <- continuous_sweep(test, fitted_svm, "skew", "raw",
                                   manual_thres = 0.25)
cat("prev_norm_trad completed \n")
timestamp()
prev_kernelr_trad <- continuous_sweep(test, fitted_svm, "kernel", "raw",
                                      manual_thres = 0.25)
cat("prev_norm_trad completed \n")
timestamp()
prev_kernelp_trad <- continuous_sweep(test, fitted_svm, "kernel", "prob",
                                      manual_thres = 0.25)
cat("prev_norm_trad completed \n")
timestamp()
prev_beta_trad <- continuous_sweep(test, fitted_svm, "beta", "prob",
                                   manual_thres = 0.25)
cat("prev_norm_trad completed \n")
timestamp()

# SLD
prev_sld <- SLD(test, fitted_svm)
cat("SLD completed \n")
timestamp()
# MS
prev_ms <- MS(test, fitted_svm)
cat("MS completed \n")
timestamp()
# DyS
prev_dys <- DyS(test, fitted_svm)
cat("DyS completed \n")
timestamp()

prev_true <- 
  readr::read_delim("dat/test_prevalences.txt", show_col_types = FALSE) %>%
  pull(`1`) 

fitted_svm$scores |>
  ggplot(aes(x = prob_scores, fill = .true_class)) +
  geom_density(alpha = 0.5)

fitted_svm$scores |>
  ggplot(aes(x = raw_scores, fill = .true_class)) +
  geom_density(alpha = 0.5)

all_prevs <-
  bind_cols(prev_norm, prev_skew, prev_kernelr, prev_kernelp, prev_beta, 
            prev_norm_trad, prev_skew_trad, prev_kernelr_trad, prev_kernelp_trad, 
            prev_beta_trad, prev_sld, prev_ms, prev_dys, prev_true)

colnames(all_prevs) <- c("cs_norm_opt", "cs_skew_opt", "cs_kernel_raw_opt", "cs_kernel_prob_opt", 
                         "cs_beta_opt", "cs_norm_trad", "cs_skew_trad", 
                         "cs_kernel_raw_trad", "cs_kernel_prob_trad", 
                         "cs_beta_trad", "sld", "ms", "dys", "prev_true")
write_rds(all_prevs, "all_prevs.RData")


# all_prevs <- readRDS("all_prevs.RData")

all_prevs |>
  pivot_longer(cs_norm_opt:dys, names_to = "quantifier", values_to = "prev_est") |>
  mutate(err = prev_est - prev_true,
         abs_err = abs(prev_est - prev_true),
         rel_err = abs_err / prev_true) |>
  ggplot(aes(x = quantifier, y = err)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

all_prevs |>
  pivot_longer(cs_norm_opt:dys, names_to = "quantifier", values_to = "prev_est") |>
  mutate(ae = abs(prev_est - prev_true),
         abs_err = abs(prev_est - prev_true),
         rel_err = abs_err / prev_true) |>
  ggplot(aes(x = quantifier, y = ae)) +  
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

compute_rel_err <- function(prev_true, prev_est, eps) {
  smooth_true <- (prev_true + eps) / (eps + 1)
  smooth_est <- (prev_est + eps) / (eps + 1)
  out1 <- abs(smooth_true - smooth_est) / smooth_true
  
  smooth_true_rev <- (1 - prev_true + eps) / (eps + 1)
  smooth_est_rev <- (1 - prev_est + eps) / (eps + 1)
  out2 <- abs(smooth_true_rev - smooth_est_rev) / smooth_true_rev
  
  out <- mean(c(out1, out2))
  return(out)
}

fitted_svm$scores |>
  ggplot(aes(x = prob_scores, fill = .true_class)) +
  geom_density(alpha = 0.6)

fitted_svm$scores |>
  ggplot(aes(x = raw_scores, fill = .true_class)) +
  geom_density(alpha = 0.6)

all_prevs |>
  pivot_longer(cs_norm_opt:dys, names_to = "quantifier", values_to = "prev_est") |>
  mutate(err = prev_est - prev_true,
         abs_err = abs(prev_est - prev_true),
         rel_err = map2_dbl(prev_true, prev_est, ~ compute_rel_err(.x, .y, 1/500))) |>
  group_by(quantifier) |>
  summarise(mae = mean(abs_err),
            rmse = sqrt(mean(abs_err^2)),
            rae = mean(rel_err)) |>
  arrange(mae)

########################################

all_prevs |>
  mutate(id = 1:n(), .before = cs_norm_opt,
         ks_dis = ks_vals) |>
  select(id, cs_skew_opt, ms, dys, sld, prev_true, ks_dis) |>
  rename(cs = cs_skew_opt) |>
  mutate(ae_cs = abs(cs - prev_true),
         ae_ms = abs(ms - prev_true),
         ae_dys = abs(dys - prev_true),
         ae_sld = abs(sld - prev_true),
         ms_cs = ae_ms - ae_cs,
         dys_cs = ae_dys - ae_cs,
         sld_cs = ae_sld - ae_cs) |>
  select(id, ks_dis, ms_cs, dys_cs, sld_cs) |>
  pivot_longer(ends_with("cs"), names_to = "quant_comp", values_to =  "diff") |>
  mutate(quant_comp = str_remove(quant_comp, "_cs")) |>
  ggplot(aes(x = ks_dis, y = diff)) +
  geom_hline(yintercept = 0) +
  geom_point(color = "grey", alpha = 0.3) +
  facet_wrap(~ quant_comp)


all_prevs |>
  mutate(id = 1:n(), .before = cs_norm_opt,
         ks_dis = ks_vals_r) |>
  select(id, cs_skew_opt, ms, dys, sld, prev_true, ks_dis) |>
  rename(cs = cs_skew_opt) |>
  mutate(ae_cs = abs(cs - prev_true),
         ae_ms = abs(ms - prev_true),
         ae_dys = abs(dys - prev_true),
         ae_sld = abs(sld - prev_true),
         ms_cs = ae_ms - ae_cs,
         dys_cs = ae_dys - ae_cs,
         sld_cs = ae_sld - ae_cs) |>
  select(id, ks_dis, ae_cs, ae_ms, ae_dys, ae_sld) |>
  pivot_longer(starts_with("ae"), names_to = "quant_comp", values_to =  "ae") |>
  mutate(quant_comp = str_remove(quant_comp, "ae_")) |>
  ggplot(aes(x = ks_dis, y = ae)) +
  geom_point(color = "grey", alpha = 0.3) +
  geom_smooth() +
  facet_wrap(~ quant_comp)


ae_dfr <- 
  all_prevs |>
  mutate(id = 1:n(), .before = cs_norm_opt,
         ks_dis = ks_vals_r) |>
  select(id, cs_skew_opt, ms, dys, sld, prev_true, ks_dis) |>
  rename(cs = cs_skew_opt) |>
  transmute(ae_cs = abs(cs - prev_true),
            ae_ms = abs(ms - prev_true),
            ae_dys = abs(dys - prev_true),
            ae_sld = abs(sld - prev_true)) |>
  pivot_longer(cols = everything(), names_to = "quantifier", values_to = "ae")

pairwise.wilcox.test(ae_dfr$ae, ae_dfr$quantifier, p.adjust.method="holm")
