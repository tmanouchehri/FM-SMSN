#############################new code

# ========================== Caporin VaR Tables ==========================
# Produces, for each dataset (Apple, S&P 500):
#  (A) Model table at 1% and 5% with columns: UC p, CC p, mean tick loss, mean viol-only loss
#  (B) Pairwise loss-differential tests (viol-only loss): upper triangle 1%, lower 5%
#      Test: z = mean(D_AB) / sqrt(Var_hat), using NeweyWest on lm(D ~ 1)
# =======================================================================

# install.packages(c("readxl","dplyr","tibble","purrr","writexl","sandwich","lmtest"))
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tibble); library(purrr)
  library(writexl); library(sandwich); library(lmtest)
})

# ---------- Helpers ----------
# Quantile (tick) loss at level tau
tick_loss <- function(y, q, tau) {
  ((y < q) * (tau - 1) + (y >= q) * tau) * (y - q)
}

# Violations-only loss (as in the note)
viol_loss <- function(y, q) {
  v <- as.integer(y < q)
  v * (1 + (y - q)^2)
}

# Kupiec unconditional coverage test
kupiec_uc <- function(viol_vec, tau) {
  n <- length(viol_vec); V <- sum(viol_vec)
  if (n == 0) return(list(rate = NA_real_, LRuc = NA_real_, p_uc = NA_real_, V = V, n = n))
  pi_hat <- V / n
  eps <- 1e-12; slog <- function(x) log(pmin(pmax(x, eps), 1 - eps))
  ll_null <- V * slog(tau) + (n - V) * slog(1 - tau)
  ll_alt  <- V * slog(pi_hat) + (n - V) * slog(1 - pi_hat)
  LRuc <- -2 * (ll_null - ll_alt)
  p_uc <- 1 - pchisq(LRuc, df = 1)
  list(rate = pi_hat, LRuc = as.numeric(LRuc), p_uc = as.numeric(p_uc), V = V, n = n)
}

# Christoffersen independence + conditional coverage
christoffersen_ind_cc <- function(viol_vec, tau) {
  v <- as.integer(viol_vec); n <- length(v)
  if (n < 2) return(list(LRind = NA_real_, p_ind = NA_real_, LRcc = NA_real_, p_cc = NA_real_))
  n00 <- n01 <- n10 <- n11 <- 0L
  for (t in 2:n) {
    if (v[t-1]==0L && v[t]==0L) n00 <- n00 + 1L
    if (v[t-1]==0L && v[t]==1L) n01 <- n01 + 1L
    if (v[t-1]==1L && v[t]==0L) n10 <- n10 + 1L
    if (v[t-1]==1L && v[t]==1L) n11 <- n11 + 1L
  }
  n0 <- n00 + n01; n1 <- n10 + n11
  p01 <- if (n0>0) n01/n0 else 0; p11 <- if (n1>0) n11/n1 else 0
  p <- if ((n0 + n1) > 0) (n01 + n11) / (n0 + n1) else 0
  eps <- 1e-12; slog <- function(x) log(pmin(pmax(x, eps), 1 - eps))
  ll_ind <- (n01 * slog(p) + n00 * slog(1 - p)) + (n11 * slog(p) + n10 * slog(1 - p))
  ll_dep <- (n01 * slog(p01) + n00 * slog(1 - p01)) + (n11 * slog(p11) + n10 * slog(1 - p11))
  LRind <- -2 * (ll_ind - ll_dep)
  p_ind <- 1 - pchisq(LRind, df = 1)
  uc <- kupiec_uc(v, tau)
  LRcc <- uc$LRuc + LRind
  p_cc <- 1 - pchisq(LRcc, df = 2)
  list(LRind = as.numeric(LRind), p_ind = as.numeric(p_ind),
       LRcc = as.numeric(LRcc), p_cc = as.numeric(p_cc))
}

# Newey–West z-test on mean(d) = 0 (regression on constant)
nw_ztest <- function(d, lag = NULL) {
  d <- as.numeric(d); d <- d[is.finite(d)]
  n <- length(d); if (n < 5) return(c(z = NA_real_, p = NA_real_))
  if (is.null(lag)) lag <- floor(1.5 * n^(1/3))
  m <- lm(d ~ 1)
  V <- sandwich::NeweyWest(m, lag = lag, prewhite = FALSE, adjust = TRUE)
  se <- sqrt(diag(V))[1]
  z  <- coef(m)[1] / se
  p  <- 2 * (1 - pnorm(abs(z)))  # asymptotic N(0,1)
  c(z = unname(z), p = unname(p))
}

# Compute VaR and all series we need
# Set is_sigma <- TRUE only if h_fore is already sigma_t (stdev). Default assumes variance.
is_sigma <- FALSE

build_var_series <- function(df) {
  names(df) <- gsub("\\s+", "", names(df))
  need <- c("omega","h_fore","y_out","Q_01","Q_05")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop(sprintf("Missing columns: %s", paste(miss, collapse=", ")))
  y <- as.numeric(df$y_out)
  s <- if (is_sigma) as.numeric(df$h_fore) else sqrt(as.numeric(df$h_fore))
  VaR1 <- as.numeric(df$omega) + s * as.numeric(df$Q_01)
  VaR5 <- as.numeric(df$omega) + s * as.numeric(df$Q_05)
  list(
    y = y, VaR1 = VaR1, VaR5 = VaR5,
    viol1 = as.integer(y < VaR1),
    viol5 = as.integer(y < VaR5),
    tick1 = tick_loss(y, VaR1, 0.01),
    tick5 = tick_loss(y, VaR5, 0.05),
    vloss1 = viol_loss(y, VaR1),
    vloss5 = viol_loss(y, VaR5)
  )
}

# Per-model metrics at a given tau
model_metrics_at_tau <- function(sers, tau = 0.01) {
  if (tau == 0.01) {
    viol <- sers$viol1; tick <- sers$tick1; vloss <- sers$vloss1
  } else {
    viol <- sers$viol5; tick <- sers$tick5; vloss <- sers$vloss5
  }
  uc <- kupiec_uc(viol, tau)
  cc <- christoffersen_ind_cc(viol, tau)
  tibble(
    V = uc$V, n = uc$n, viol_rate = uc$rate,
    UC_p = uc$p_uc, CC_p = cc$p_cc,
    mean_tick_loss = mean(tick,  na.rm=TRUE),
    mean_viol_loss = mean(vloss, na.rm=TRUE)
  )
}

# Build the requested tables for one dataset
make_caporin_tables <- function(files, dataset_label = "") {
  # Read and build series
  dfs <- lapply(files, read_excel)
  names(dfs) <- basename(files)
  ser <- lapply(dfs, build_var_series)
  
  # (A) Model tables at 1% and 5%
  tab1 <- imap_dfr(ser, ~{
    m <- model_metrics_at_tau(.x, 0.01)
    tibble(dataset = dataset_label, model = .y, m)
  }) %>% select(dataset, model, UC_p, CC_p, mean_tick_loss, mean_viol_loss, viol_rate, V, n)
  
  tab5 <- imap_dfr(ser, ~{
    m <- model_metrics_at_tau(.x, 0.05)
    tibble(dataset = dataset_label, model = .y, m)
  }) %>% select(dataset, model, UC_p, CC_p, mean_tick_loss, mean_viol_loss, viol_rate, V, n)
  
  # (B) Pairwise z-tests for violations-only loss differentials
  mods <- names(ser); M <- length(mods)
  Zmat <- matrix(NA_real_, M, M, dimnames = list(mods, mods))
  Pmat <- matrix(NA_real_, M, M, dimnames = list(mods, mods))
  
  for (i in seq_len(M)) {
    for (j in seq_len(M)) {
      if (i == j) next
      if (i < j) {
        # Upper triangle: tau = 1%
        d <- ser[[mods[i]]]$vloss1 - ser[[mods[j]]]$vloss1
        st <- nw_ztest(d)
        Zmat[i, j] <- st["z"]; Pmat[i, j] <- st["p"]
      } else {
        # Lower triangle: tau = 5%
        d <- ser[[mods[i]]]$vloss5 - ser[[mods[j]]]$vloss5
        st <- nw_ztest(d)
        Zmat[i, j] <- st["z"]; Pmat[i, j] <- st["p"]
      }
    }
  }
  
  list(table_1pct = tab1, table_5pct = tab5, pairwise_Z = Zmat, pairwise_p = Pmat)
}

# ===================== RUN FOR YOUR FILES ======================
apple_files <- c(
  "normal_garch_results_apple.xlsx",
  "mixture_garch_results_apple.xlsx",
  "mixture_garch_results_st_applet.xlsx"
)

sp_files <- c(
  "normal_garch_results.xlsx",
  "mixture_garch_results.xlsx",
  "mixture_garch_results_st.xlsx"
)

apple_out <- make_caporin_tables(apple_files, "Apple")
sp_out    <- make_caporin_tables(sp_files, "S&P 500")

# Inspect in console (examples)
apple_out$table_1pct
apple_out$table_5pct
apple_out$pairwise_Z   # z-statistics: upper=1% , lower=5%
apple_out$pairwise_p   # p-values:    upper=1% , lower=5%

sp_out$table_1pct
sp_out$table_5pct
sp_out$pairwise_Z
sp_out$pairwise_p

# ===================== WRITE TO EXCEL ==========================
write_xlsx(
  list(
    Apple_Table_1pct = apple_out$table_1pct,
    Apple_Table_5pct = apple_out$table_5pct,
    Apple_Pairwise_Z = as.data.frame(apple_out$pairwise_Z),
    Apple_Pairwise_p = as.data.frame(apple_out$pairwise_p),
    
    SP500_Table_1pct = sp_out$table_1pct,
    SP500_Table_5pct = sp_out$table_5pct,
    SP500_Pairwise_Z = as.data.frame(sp_out$pairwise_Z),
    SP500_Pairwise_p = as.data.frame(sp_out$pairwise_p)
  ),
  path = "VaR_Caporin_Tables.xlsx"
)

# ===================== INTERPRETATION NOTE =====================
# Pairwise table entries are z-statistics for mean(L_A - L_B) with NW SE.
# If z < 0 and the p-value is small (say < 0.05), model A has LOWER loss than B
# (i.e., A is preferred). Upper triangle uses τ=1%; lower triangle uses τ=5%.




# One-page Apple VaR comparison (colored 1% vs 5%)
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
})

# If h_fore is already sigma_t (stdev) set TRUE; otherwise it's variance (default)
is_sigma <- FALSE

compute_var_tidy <- function(path, model_label) {
  d <- readxl::read_excel(path)
  names(d) <- gsub("\\s+", "", names(d))
  need <- c("omega","h_fore","y_out","Q_01","Q_05")
  miss <- setdiff(need, names(d)); if (length(miss)) stop(paste("Missing:", paste(miss, collapse=", ")))
  
  # use date if present else index
  date_cols <- grep("date|time|idx|index", names(d), ignore.case = TRUE, value = TRUE)
  t <- if (length(date_cols)) d[[date_cols[1]]] else seq_len(nrow(d))
  
  s <- if (is_sigma) as.numeric(d$h_fore) else sqrt(as.numeric(d$h_fore))
  VaR01 <- as.numeric(d$omega) + s * as.numeric(d$Q_01)
  VaR05 <- as.numeric(d$omega) + s * as.numeric(d$Q_05)
  y     <- as.numeric(d$y_out)
  
  tibble(
    t = t,
    Return = y,
    `VaR 1%` = VaR01,
    `VaR 5%` = VaR05,
    Model = model_label
  ) |>
    pivot_longer(c(`VaR 1%`,`VaR 5%`), names_to = "VaR_level", values_to = "VaR") |>
    mutate(Violation = Return < VaR)
}

# File -> label
apple_files <- c(
  "normal_garch_results_apple.xlsx"       = "N–GARCH",
  "mixture_garch_results_apple.xlsx"      = "FM2–N–GARCH",
  "mixture_garch_results_st_applet.xlsx"  = "FM2–ST–GARCH"
)

apple_all <- bind_rows(lapply(names(apple_files), function(p)
  compute_var_tidy(p, apple_files[[p]])
))

p <- ggplot(apple_all, aes(x = t)) +
  geom_line(aes(y = Return), linewidth = 0.4, color = "grey35") +
  geom_line(aes(y = VaR, color = VaR_level, linetype = VaR_level), linewidth = 0.8) +
  geom_point(
    data = subset(apple_all, Violation),
    aes(y = Return, color = VaR_level),
    size = 1.5, alpha = 0.9
  ) +
  scale_color_manual(
    values = c("VaR 1%" = "#D62728",  # red
               "VaR 5%" = "#1F77B4"), # blue
    name = "VaR level"
  ) +
  scale_linetype_manual(
    values = c("VaR 1%" = "dashed", "VaR 5%" = "dotdash"),
    name = "VaR level"
  ) +
  labs(title = "Apple — VaR (1% & 5%) across models",
       x = NULL, y = "Return / VaR") +
  facet_wrap(~ Model, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"))

# Show in R
print(p)

# Save one-page figure (PDF and PNG)
ggsave("Apple_VaR_all_models_onepage.pdf", p, width = 9, height = 10)
ggsave("Apple_VaR_all_models_onepage.png", p, width = 9, height = 10, dpi = 220)
