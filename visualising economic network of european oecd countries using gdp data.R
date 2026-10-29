############################################################
## European OECD GDP Time-Varying Network Dashboard
## Methodology: factor-augmented endogenous network model
## with CCE control-function construction
##
## Outcome:
##   Y_it = ln_gdp_percap
##
## Auxiliary network variable:
##   Z_it = ln_openness
##
## Auxiliary equation estimated by CCE-OLS:
##   Z_it = beta_iz' X_itz + sum_{l=0}^{pT} pi_i,lz' rbar_{t-l} + u_it
##
## Outcome equation estimated by CFQR-CCE:
##   Y_it = alpha_i + Gamma_i [W_t Y_t]_i + lambda_i Y_i,t-1
##          + rho_i [W_{t-1}Y_{t-1}]_i + beta_iy'X_ity
##          + sum_{l=0}^{pT} Delta_i,l' rbar_{t-l}
##          + delta_i epshat_it + error_it
##
## Network:
##   W_ij,t proportional to exp{-eta |Z_it - Z_jt|} 1(d_ij <= d_c),
##   where d_ij is computed from country longitude/latitude coordinates.
############################################################

rm(list = ls())

############################################################
## Required packages
############################################################

required_packages <- c(
  "shiny", "shinydashboard", "plotly", "dplyr", "tidyr", "readr", "scales", "quantreg"
)

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Please install the following packages before running this app: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(shiny)
library(shinydashboard)
library(plotly)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(quantreg)

#########################################################################
## Selct and assign the variable for Y_{it}, \bfX_{ity}, Z_{it}, \bfX_{itz}
###########################################################################

DATA_FILE <- "GDP data.csv"
TAU <- 0.50
ETA_WEIGHT <- 1.00
WEIGHT_THRESHOLD <- 0.02
DC_DEGREES <- 15
PT_DEFAULT <- 2

Y_VAR <- "ln_gdp_percap"
Z_VAR <- "ln_openness" # log [(imports + exports)/ GDP]

XZ_VARS <- c(
  "fdi", # log (foreign direct investment)
  "ln_imports", 
  "ln_exports",
  "ln_rnd_imp_spill"
)

XY_VARS <- c(
  "ln_physical_capital",
  "ln_labour_force",
  "ln_human_capital",
  "ln_domestic_rnd",
  "ln_infrastructure",
  "ln_urbanisation",
  "pol_stability",
  "corruption"
)

MODEL_VARS <- c(Y_VAR, Z_VAR, XZ_VARS, XY_VARS)

######################################################################
## OECD member countries and longitude/latitude coordinates
## Coordinates are approximated by considering the country coordinates of the capitals
# for each country
########################################################################

oecd_members <- c(
  "Australia", "Austria", "Belgium", "Canada", "Chile", "Colombia",
  "Costa Rica", "Czech Republic", "Denmark", "Estonia", "Finland",
  "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland",
  "Israel", "Italy", "Japan", "Republic of Korea", "Latvia", "Lithuania",
  "Luxembourg", "Mexico", "Netherlands", "New Zealand", "Norway",
  "Poland", "Portugal", "Slovak Republic", "Slovenia", "Spain",
  "Sweden", "Switzerland", "Turkiye", "United Kingdom", "United States"
)

european_oecd_members <- c(
  "Austria", "Belgium", "Czech Republic", "Denmark", "Estonia", "Finland",
  "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy",
  "Latvia", "Lithuania", "Luxembourg", "Netherlands", "Norway", "Poland",
  "Portugal", "Slovak Republic", "Slovenia", "Spain", "Sweden", "Switzerland",
  "Turkiye", "United Kingdom"
)

country_coords <- data.frame(
  country = c(
    "Australia", "Austria", "Belgium", "Canada", "Chile", "Colombia",
    "Costa Rica", "Czech Republic", "Denmark", "Estonia", "Finland",
    "France", "Germany", "Greece", "Hungary", "Iceland", "Ireland",
    "Israel", "Italy", "Japan", "Republic of Korea", "Latvia", "Lithuania",
    "Luxembourg", "Mexico", "Netherlands", "New Zealand", "Norway",
    "Poland", "Portugal", "Slovak Republic", "Slovenia", "Spain",
    "Sweden", "Switzerland", "Turkiye", "United Kingdom", "United States"
  ),
  longitude = c(
    133.7751, 14.5501, 4.4699, -106.3468, -71.5430, -74.2973,
    -83.7534, 15.4730, 9.5018, 25.0136, 25.7482,
    2.2137, 10.4515, 21.8243, 19.5033, -19.0208, -8.2439,
    34.8516, 12.5674, 138.2529, 127.7669, 24.6032, 23.8813,
    6.1296, -102.5528, 5.2913, 174.8860, 8.4689,
    19.1451, -8.2245, 19.6990, 14.9955, -3.7492,
    18.6435, 8.2275, 35.2433, -3.4360, -95.7129
  ),
  latitude = c(
    -25.2744, 47.5162, 50.5039, 56.1304, -35.6751, 4.5709,
    9.7489, 49.8175, 56.2639, 58.5953, 61.9241,
    46.2276, 51.1657, 39.0742, 47.1625, 64.9631, 53.4129,
    31.0461, 41.8719, 36.2048, 35.9078, 56.8796, 55.1694,
    49.8153, 23.6345, 52.1326, -40.9006, 60.4720,
    51.9194, 39.3999, 48.6690, 46.1512, 40.4637,
    60.1282, 46.8182, 38.9637, 55.3781, 37.0902
  ),
  stringsAsFactors = FALSE
)

# Defensive check: dplyr::left_join(..., by = "country") needs this exact column.
stopifnot(all(c("country", "longitude", "latitude") %in% names(country_coords)))

############################################################
## scaling to to draw the nodes and the edges and lags for the dynamic terms
############################################################

safe_scale_vec <- function(x) {
  x <- as.numeric(x)
  sx <- sd(x, na.rm = TRUE)
  mx <- mean(x, na.rm = TRUE)
  if (!is.finite(sx) || sx < 1e-12) return(rep(0, length(x)))
  (x - mx) / sx
}

scale_range <- function(x, a, b) {
  x <- as.numeric(x)
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  if (!is.finite(xmin) || !is.finite(xmax) || abs(xmax - xmin) < 1e-12) {
    return(rep((a + b) / 2, length(x)))
  }
  a + (x - xmin) / (xmax - xmin) * (b - a)
}

make_lag <- function(x, lag = 1) {
  if (lag <= 0) return(x)
  c(rep(NA_real_, lag), x[1:(length(x) - lag)])
}

############################################################
## Distance matrix from longitude/latitude coordinates
############################################################

make_geo_linf_distance_matrix <- function(lon, lat) {
  N <- length(lon)
  D <- matrix(0, N, N)
  for (i in 1:N) {
    for (j in 1:N) {
      dlon_raw <- abs(lon[i] - lon[j])
      dlon <- min(dlon_raw, 360 - dlon_raw)
      dlat <- abs(lat[i] - lat[j])
      D[i, j] <- max(dlon, dlat)
    }
  }
  D
}

############################################################
## CCE nodal averages: rbar_t = (Ybar_t, Xbar_y,t', Zbar_t, Xbar_z,t')'
############################################################

make_rbar <- function(panel_df, time_values, y_var, z_var, xy_vars, xz_vars) {
  rbar <- panel_df %>%
    group_by(year_actual) %>%
    summarise(
      Ybar = mean(.data[[y_var]], na.rm = TRUE),
      across(all_of(xy_vars), ~mean(.x, na.rm = TRUE), .names = "Yavg_{.col}"),
      Zbar = mean(.data[[z_var]], na.rm = TRUE),
      across(all_of(xz_vars), ~mean(.x, na.rm = TRUE), .names = "Zavg_{.col}"),
      .groups = "drop"
    ) %>%
    arrange(year_actual)

  rbar_mat <- as.matrix(rbar %>% select(-year_actual))
  rownames(rbar_mat) <- rbar$year_actual
  rbar_mat
}

make_rbar_lag_matrix <- function(rbar_mat, pT) {
  Tn <- nrow(rbar_mat)
  out <- NULL
  out_names <- character(0)
  for (ell in 0:pT) {
    lagged <- matrix(NA_real_, Tn, ncol(rbar_mat))
    if (ell == 0) {
      lagged <- rbar_mat
    } else {
      lagged[(ell + 1):Tn, ] <- rbar_mat[1:(Tn - ell), , drop = FALSE]
    }
    colnames(lagged) <- paste0(colnames(rbar_mat), "_L", ell)
    out <- cbind(out, lagged)
    out_names <- c(out_names, colnames(lagged))
  }
  colnames(out) <- out_names
  rownames(out) <- rownames(rbar_mat)
  out
}

############################################################
## Estimate Z equation by unit-specific CCE-OLS
############################################################

estimate_z_cce_ols <- function(panel_df, countries, time_values, z_var, xz_vars, rbar_lags) {
  out <- vector("list", length(countries))
  names(out) <- countries

  for (cc in countries) {
    df_i <- panel_df %>% filter(country == cc) %>% arrange(year_actual)
    R_i <- rbar_lags[as.character(df_i$year_actual), , drop = FALSE]

    est_df <- bind_cols(
      df_i %>% select(year_actual, all_of(z_var), all_of(xz_vars)),
      as.data.frame(R_i)
    ) %>%
      tidyr::drop_na()

    if (nrow(est_df) <= length(xz_vars) + 3) {
      df_i$epshat <- NA_real_
      out[[cc]] <- list(data = df_i, fit = NULL, coef = NA)
      next
    }

    form <- as.formula(
      paste(z_var, "~", paste(c(xz_vars, colnames(R_i)), collapse = " + "))
    )

    fit <- tryCatch(lm(form, data = est_df), error = function(e) NULL)

    df_i$epshat <- NA_real_
    if (!is.null(fit)) {
      est_df$epshat <- residuals(fit)
      df_i <- df_i %>%
        left_join(est_df %>% select(year_actual, epshat), by = "year_actual") %>%
        mutate(epshat = coalesce(epshat.y, epshat.x)) %>%
        select(-epshat.x, -epshat.y)
    }

    out[[cc]] <- list(data = df_i, fit = fit, coef = if (!is.null(fit)) coef(fit) else NA)
  }

  bind_rows(lapply(out, `[[`, "data"))
}

############################################################
## Build W_t from observed Z_it and geographic distance
############################################################

# h(Z_{it}, Z_{jt}) = exp(-1*|Z_{it}- Z_{jt}|) example from Baltagi and Shu 2025
h_kernel <- function(z_i, z_j, eta = 1.0) { 
  exp(-eta * abs(z_i - z_j))
}

build_Wt_geo <- function(
    z_t,
    Dmat,
    dc = DC_DEGREES,
    eta = ETA_WEIGHT,
    weight_threshold = WEIGHT_THRESHOLD,
    renormalize_after_threshold = TRUE
) {
  N <- length(z_t)

  H <- outer(z_t, z_t, Vectorize(function(a, b) h_kernel(a, b, eta = eta)))

  local_mask <- 1 * (Dmat <= dc)
  diag(local_mask) <- 0

  Wraw <- H * local_mask
  diag(Wraw) <- 0

  rs <- rowSums(Wraw)
  W <- matrix(0, N, N)
  idx <- rs > 0
  W[idx, ] <- Wraw[idx, , drop = FALSE] / rs[idx]

  W[W < weight_threshold] <- 0
  diag(W) <- 0

  if (renormalize_after_threshold) {
    rs2 <- rowSums(W)
    idx2 <- rs2 > 0
    W[idx2, ] <- W[idx2, , drop = FALSE] / rs2[idx2]
  }

  W
}

build_Wlist_geo <- function(Z_mat, Dmat, dc = DC_DEGREES, eta = ETA_WEIGHT, weight_threshold = WEIGHT_THRESHOLD) {
  Tn <- ncol(Z_mat)
  Wlist <- vector("list", Tn)
  for (tt in 1:Tn) {
    Wlist[[tt]] <- build_Wt_geo(
      z_t = Z_mat[, tt],
      Dmat = Dmat,
      dc = dc,
      eta = eta,
      weight_threshold = weight_threshold,
      renormalize_after_threshold = TRUE
    )
  }
  Wlist
}

############################################################
## Estimate Y equation by unit-specific CFQR-CCE
############################################################

estimate_y_cfqr_cce <- function(panel_df, countries, time_values, y_var, xy_vars, rbar_lags, Wlist, tau = TAU) {
  fit_list <- vector("list", length(countries))
  names(fit_list) <- countries

  Y_mat <- xtabs(as.formula(paste(y_var, "~ country + year_actual")), data = panel_df)
  country_order <- rownames(Y_mat)
  year_order <- as.integer(colnames(Y_mat))

  for (idx_country in seq_along(countries)) {
    cc <- countries[idx_country]
    i <- match(cc, country_order)
    if (is.na(i)) next

    df_i <- panel_df %>% filter(country == cc) %>% arrange(year_actual)
    t_idx <- match(df_i$year_actual, year_order)

    WY <- rep(NA_real_, nrow(df_i))
    WlagYlag <- rep(NA_real_, nrow(df_i))
    Ylag <- make_lag(df_i[[y_var]], 1)

    for (r in seq_along(t_idx)) {
      tt <- t_idx[r]
      if (!is.na(tt)) {
        WY[r] <- as.numeric(Wlist[[tt]][i, ] %*% Y_mat[, tt])
      }
      if (!is.na(tt) && tt > 1) {
        WlagYlag[r] <- as.numeric(Wlist[[tt - 1]][i, ] %*% Y_mat[, tt - 1])
      }
    }

    R_i <- rbar_lags[as.character(df_i$year_actual), , drop = FALSE]

    est_df <- bind_cols(
      df_i %>% select(year_actual, all_of(y_var), all_of(xy_vars), epshat),
      tibble(WY = WY, Ylag = Ylag, WlagYlag = WlagYlag),
      as.data.frame(R_i)
    ) %>%
      tidyr::drop_na()

    if (nrow(est_df) <= length(xy_vars) + 6) {
      # Keep the list element and the country name. Do not use [[cc]] <- NULL,
      # because that deletes the element and later removes the country column
      # from the coefficient summary.
      fit_list[idx_country] <- list(NULL)
      next
    }

    form <- as.formula(
      paste(
        y_var, "~",
        paste(c("WY", "Ylag", "WlagYlag", xy_vars, colnames(R_i), "epshat"), collapse = " + ")
      )
    )

    fit <- tryCatch(
      quantreg::rq(form, tau = tau, data = est_df, method = "br"),
      error = function(e) NULL
    )

    # Keep the list element even when the fitted object is NULL.
    fit_list[idx_country] <- list(fit)
  }

  names(fit_list) <- countries
  fit_list
}

extract_main_coefficients <- function(fit_list, countries) {
  # Always return a tibble with a country column, even if all rq fits fail.
  safe_coef <- function(cf, nm) {
    if (is.null(cf) || !(nm %in% names(cf)) || is.na(cf[[nm]])) {
      return(NA_real_)
    }
    as.numeric(cf[[nm]])
  }

  rows <- lapply(countries, function(cc) {
    fit <- fit_list[[cc]]

    if (is.null(fit)) {
      return(tibble(
        country = as.character(cc),
        Gamma_WY = NA_real_,
        lambda_Ylag = NA_real_,
        rho_WlagYlag = NA_real_,
        delta_epshat = NA_real_
      ))
    }

    cf <- stats::coef(fit)

    tibble(
      country = as.character(cc),
      Gamma_WY = safe_coef(cf, "WY"),
      lambda_Ylag = safe_coef(cf, "Ylag"),
      rho_WlagYlag = safe_coef(cf, "WlagYlag"),
      delta_epshat = safe_coef(cf, "epshat")
    )
  })

  out <- bind_rows(rows)

  # Defensive check: this prevents the left_join error seen previously.
  if (!"country" %in% names(out)) {
    out <- tibble(
      country = as.character(countries),
      Gamma_WY = NA_real_,
      lambda_Ylag = NA_real_,
      rho_WlagYlag = NA_real_,
      delta_epshat = NA_real_
    )
  }

  out
}

############################################################
## Load and prepare OECD data
############################################################

if (!file.exists(DATA_FILE)) {
  stop("Could not find GDP data.csv in the app folder. Please place GDP data.csv in the same folder as this app.R file.", call. = FALSE)
}

raw <- readr::read_csv(DATA_FILE, show_col_types = FALSE) %>%
  mutate(country = trimws(as.character(country)))

country_coords <- country_coords %>%
  mutate(country = trimws(as.character(country)))

needed_cols <- c("country", "year_actual", MODEL_VARS)
missing_cols <- setdiff(needed_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("The dataset is missing the following required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

panel0 <- raw %>%
  filter(country %in% european_oecd_members) %>%
  dplyr::left_join(country_coords, by = "country") %>%
  filter(!is.na(longitude), !is.na(latitude)) %>%
  select(country, iso_code, region, year_actual, longitude, latitude, all_of(MODEL_VARS)) %>%
  arrange(country, year_actual)

complete_countries <- panel0 %>%
  group_by(country) %>%
  summarise(
    n_years = n_distinct(year_actual),
    min_year = min(year_actual),
    max_year = max(year_actual),
    complete_vars = all(stats::complete.cases(across(all_of(MODEL_VARS)))),
    .groups = "drop"
  ) %>%
  filter(complete_vars) %>%
  pull(country)

common_years <- panel0 %>%
  filter(country %in% complete_countries) %>%
  count(year_actual) %>%
  filter(n == length(complete_countries)) %>%
  arrange(year_actual) %>%
  pull(year_actual)

panel <- panel0 %>%
  filter(country %in% complete_countries, year_actual %in% common_years) %>%
  arrange(country, year_actual)

countries <- sort(unique(panel$country))
years <- sort(unique(panel$year_actual))
N <- length(countries)
Tn <- length(years)

if (N < 3 || Tn < 5) {
  stop("After filtering to OECD countries and complete cases, too few observations remain.", call. = FALSE)
}

pT <- min(PT_DEFAULT, max(0, floor((Tn - 5) / 4)))

############################################################
## Standardized variables for numerical stability for CCE
## we use original values for labels; use standardized variables in estimation.
############################################################

panel_std <- panel %>%
  mutate(across(all_of(MODEL_VARS), safe_scale_vec, .names = "std_{.col}"))

Y_STD <- paste0("std_", Y_VAR)
Z_STD <- paste0("std_", Z_VAR)
XZ_STD <- paste0("std_", XZ_VARS)
XY_STD <- paste0("std_", XY_VARS)

############################################################
## CCE first stage: estimate Z equation and residual epshat
############################################################

rbar_mat <- make_rbar(
  panel_df = panel_std,
  time_values = years,
  y_var = Y_STD,
  z_var = Z_STD,
  xy_vars = XY_STD,
  xz_vars = XZ_STD
)

rbar_lags <- make_rbar_lag_matrix(rbar_mat, pT = pT)

panel_cce <- estimate_z_cce_ols(
  panel_df = panel_std,
  countries = countries,
  time_values = years,
  z_var = Z_STD,
  xz_vars = XZ_STD,
  rbar_lags = rbar_lags
)

panel_cce <- panel_cce %>%
  group_by(country) %>%
  mutate(epshat = ifelse(is.na(epshat), mean(epshat, na.rm = TRUE), epshat)) %>%
  ungroup() %>%
  mutate(epshat = ifelse(is.na(epshat), 0, epshat))

############################################################
## Matrices ordered country x time
############################################################

panel_cce$country <- factor(panel_cce$country, levels = countries)
panel_cce$year_actual <- factor(panel_cce$year_actual, levels = years)

Y_mat <- xtabs(as.formula(paste(Y_STD, "~ country + year_actual")), data = panel_cce)
Z_mat <- xtabs(as.formula(paste(Z_STD, "~ country + year_actual")), data = panel_cce)
Y_original_mat <- xtabs(as.formula(paste(Y_VAR, "~ country + year_actual")), data = panel_cce)
Z_original_mat <- xtabs(as.formula(paste(Z_VAR, "~ country + year_actual")), data = panel_cce)

coords_df <- panel_cce %>%
  distinct(country, longitude, latitude) %>%
  arrange(match(country, countries))

coords <- as.matrix(coords_df[, c("longitude", "latitude")])
Dmat <- make_geo_linf_distance_matrix(coords[, 1], coords[, 2])

Wlist <- build_Wlist_geo(
  Z_mat = Z_mat,
  Dmat = Dmat,
  dc = DC_DEGREES,
  eta = ETA_WEIGHT,
  weight_threshold = WEIGHT_THRESHOLD
)

############################################################
## CFQR-CCE outcome estimation
############################################################

fit_y_list <- estimate_y_cfqr_cce(
  panel_df = panel_cce,
  countries = countries,
  time_values = years,
  y_var = Y_STD,
  xy_vars = XY_STD,
  rbar_lags = rbar_lags,
  Wlist = Wlist,
  tau = TAU
)

coef_summary <- extract_main_coefficients(fit_y_list, countries)

############################################################
## Summaries for dashboard
############################################################

network_summary <- tibble(
  year_actual = years,
  edges = sapply(Wlist, function(W) sum(W > 0)),
  density = sapply(Wlist, function(W) sum(W > 0) / (N * (N - 1))),
  avg_weight = sapply(Wlist, function(W) ifelse(sum(W > 0) == 0, NA_real_, mean(W[W > 0]))),
  max_weight = sapply(Wlist, function(W) max(W)),
  avg_Y = as.numeric(colMeans(Y_original_mat)),
  sd_Y = apply(Y_original_mat, 2, sd),
  avg_Z = as.numeric(colMeans(Z_original_mat)),
  avg_epshat = panel_cce %>%
    mutate(year_actual_num = as.integer(as.character(year_actual))) %>%
    group_by(year_actual_num) %>%
    summarise(avg_epshat = mean(epshat, na.rm = TRUE), .groups = "drop") %>%
    pull(avg_epshat)
)

# Make sure both join keys are character variables.
coef_summary <- coef_summary %>% mutate(country = as.character(country))

panel_plot <- panel_cce %>%
  mutate(
    country = as.character(country),
    year_actual = as.integer(as.character(year_actual)),
    Y = .data[[Y_VAR]],
    Z = .data[[Z_VAR]],
    Y_std = .data[[Y_STD]],
    Z_std = .data[[Z_STD]]
  ) %>%
  left_join(coef_summary, by = "country")

############################################################
## Plotly network visualization
############################################################

make_network_plotly <- function(W, coords, Y_t, Z_t, epshat_t, year_t, country_names, dc = DC_DEGREES) {
  N <- nrow(W)

  node_df <- tibble(
    country = country_names,
    lon = coords[, 1],
    lat = coords[, 2],
    Y = as.numeric(Y_t),
    Z = as.numeric(Z_t),
    epshat = as.numeric(epshat_t),
    node_size = scale_range(Y_t, 14, 34)
  )

  edge_idx <- which(W > 0, arr.ind = TRUE)

  # Zoom to the European OECD coordinate window so that links are visible.
  x_rng <- range(coords[, 1], na.rm = TRUE)
  y_rng <- range(coords[, 2], na.rm = TRUE)
  x_pad <- max(2.0, 0.08 * diff(x_rng))
  y_pad <- max(1.5, 0.12 * diff(y_rng))
  x_plot_range <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)
  y_plot_range <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)

  p <- plot_ly()

  if (nrow(edge_idx) > 0) {
    edge_weights <- W[edge_idx]
    edge_widths <- scale_range(edge_weights, 0.4, 3.2)

    for (k in seq_len(nrow(edge_idx))) {
      i <- edge_idx[k, 1]
      j <- edge_idx[k, 2]

      p <- p %>%
        add_trace(
          x = c(coords[j, 1], coords[i, 1], NA),
          y = c(coords[j, 2], coords[i, 2], NA),
          type = "scatter",
          mode = "lines",
          line = list(width = edge_widths[k], color = "rgba(70,70,70,0.22)"),
          hoverinfo = "text",
          text = paste0(
            "Directed link: ", country_names[j], " → ", country_names[i],
            "<br>W<sub>", i, j, ",", year_t, "</sub> = ", round(W[i, j], 3),
            "<br>Distance d<sub>ij</sub> = ", round(Dmat[i, j], 2), " degrees",
            "<br>Cutoff d<sub>c</sub> = ", dc, " degrees"
          ),
          showlegend = FALSE
        )
    }
  }

  p %>%
    add_trace(
      data = node_df,
      x = ~lon,
      y = ~lat,
      type = "scatter",
      mode = "markers+text",
      text = ~country,
      textposition = "top center",
      textfont = list(size = 9, color = "black"),
      marker = list(
        size = ~node_size,
        color = ~Y,
        colorscale = list(
          c(0, "blue"), c(0.2, "cyan"), c(0.4, "green"),
          c(0.6, "yellow"), c(0.8, "orange"), c(1, "red")
        ),
        showscale = TRUE,
        colorbar = list(title = list(text = "ln GDP<br>per capita"), thickness = 14, len = 0.75, x = 1.03),
        line = list(color = "black", width = 1)
      ),
      hovertemplate = paste0(
        "<b>%{text}</b>",
        "<br>Longitude: %{x:.2f}",
        "<br>Latitude: %{y:.2f}",
        "<br>ln GDP per capita: %{marker.color:.2f}",
        "<extra></extra>"
      ),
      showlegend = FALSE
    ) %>%
    layout(
      title = list(text = paste0("European OECD GDP Network, ", year_t), x = 0.02),
      xaxis = list(
        title = "Longitude",
        range = x_plot_range,
        zeroline = FALSE,
        gridcolor = "rgba(220,220,220,0.6)"
      ),
      yaxis = list(
        title = "Latitude",
        range = y_plot_range,
        zeroline = FALSE,
        gridcolor = "rgba(220,220,220,0.6)"
      ),
      plot_bgcolor = "rgba(255,255,255,0.95)",
      paper_bgcolor = "white",
      margin = list(l = 60, r = 120, t = 60, b = 60),
      font = list(family = "Gill Sans, Arial, sans-serif")
    ) %>%
    config(displayModeBar = FALSE)
}

############################################################
## User interface for the dashboard
############################################################

ui <- dashboardPage(
  title = "European OECD GDP Network",
  skin = "black",

  dashboardHeader(title = "European OECD GDP Network Dashboard", titleWidth = 330),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      id = "tabs",
      menuItem("Current Network", tabName = "current", icon = icon("project-diagram")),
      menuItem("Evolution", tabName = "evolution", icon = icon("chart-line")),
      menuItem("CCE/CFQR Analysis", tabName = "analysis", icon = icon("calculator"))
    ),

    div(
      style = "padding: 15px; color: #d8d8d8;",
      p("European OECD country network using longitude/latitude topology and time-varying weights from observed openness."),
      p("Node colour and size represent ln GDP per capita. Edge width represents W_ij,t.")
    ),

    div(
      style = "padding: 15px;",
      h5("Select Year"),
      sliderInput(
        "year", NULL,
        min = min(years), max = max(years), value = min(years),
        step = 1, sep = "", animate = FALSE, width = "100%", ticks = TRUE
      ),
      actionButton("playBtn", label = "", icon = icon("play"), class = "btn btn-primary")
    ),

    div(
      style = "padding: 15px;",
      h5("Network Statistics"),
      strong("Average ln GDP per capita"),
      div(style = "font-size: 25px; color: #00c0ef;", textOutput("avgYText")),
      br(),
      strong("Nonzero directed links"),
      div(style = "font-size: 25px; color: #f39c12;", textOutput("edgeCountText")),
      br(),
      strong("Countries / CCE lag order"),
      div(style = "font-size: 18px; color: #d8d8d8;", textOutput("modelSizeText"))
    ),

    div(
      style = "position:absolute; bottom:20px; padding: 15px; color:#b8b8b8;",
      p("Data: GDP data.csv; European OECD subset"),
      p("Model: CCE control-function network quantile model")
    )
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        * { font-family: 'Gill Sans', Arial, sans-serif !important; }
        body, .content-wrapper { background-color: #f7f7f7; }
        .section-header { background: white; padding: 18px 22px; margin-bottom: 15px; border-radius: 8px; box-shadow: 0 1px 5px rgba(0,0,0,0.08); }
        .section-subtitle { color: #666; font-size: 15px; margin-top: 5px; }
        .viz-container { background: white; padding: 15px; border-radius: 8px; box-shadow: 0 1px 5px rgba(0,0,0,0.08); }
        .stat-box { background: white; text-align: center; padding: 18px; border-radius: 8px; box-shadow: 0 1px 5px rgba(0,0,0,0.08); }
        .stat-value { font-size: 30px; font-weight: bold; color: #2c3e50; }
        .stat-label { font-size: 13px; color: #777; margin-top: 5px; }
        .info-box-custom { background: white; padding: 15px; border-left: 4px solid #3c8dbc; border-radius: 6px; box-shadow: 0 1px 5px rgba(0,0,0,0.08); margin-bottom: 15px; }
      ")),
      tags$script(HTML("
        $(document).ready(function() {
          var playing = false;
          var interval;
          $(document).on('click', '#playBtn', function(e) {
            e.preventDefault();
            if (!playing) {
              playing = true;
              $(this).find('i').removeClass('fa-play').addClass('fa-pause');
              var slider = $('#year').data('ionRangeSlider');
              if (!slider) return;
              var currentVal = slider.result.from;
              var maxVal = slider.result.max;
              var minVal = slider.result.min;
              if (currentVal >= maxVal) slider.update({from: minVal});
              interval = setInterval(function() {
                var current = slider.result.from;
                if (current < maxVal) {
                  slider.update({from: current + 1});
                } else {
                  clearInterval(interval);
                  playing = false;
                  $('#playBtn').find('i').removeClass('fa-pause').addClass('fa-play');
                }
              }, 1300);
            } else {
              playing = false;
              clearInterval(interval);
              $(this).find('i').removeClass('fa-pause').addClass('fa-play');
            }
          });
        });
      "))
    ),

    tabItems(
      tabItem(
        tabName = "current",
        fluidRow(
          column(12, div(class = "section-header",
            h3("Current European OECD GDP Network"),
            div(class = "section-subtitle", "Longitude/latitude topology with weights determined by openness similarity"),
            p("Each node is an OECD country. The outcome ", tags$b("Y_it"), " is ", tags$b("ln_gdp_percap"), ". The auxiliary variable ", tags$b("Z_it"), " is ", tags$b("ln_openness"), ".")
          ))
        ),
        fluidRow(column(12, div(class = "viz-container", plotlyOutput("networkPlot", height = "650px")))),
        fluidRow(style = "margin-top: 15px;",
          column(3, div(class = "stat-box", div(class = "stat-value", textOutput("statCountries")), div(class = "stat-label", "OECD Countries in Dataset"))),
          column(3, div(class = "stat-box", div(class = "stat-value", textOutput("statEdges")), div(class = "stat-label", "Nonzero Directed Links"))),
          column(3, div(class = "stat-box", div(class = "stat-value", textOutput("statDensity")), div(class = "stat-label", "Network Density"))),
          column(3, div(class = "stat-box", div(class = "stat-value", textOutput("statAvgOutcome")), div(class = "stat-label", "Average ln GDP per capita")))
        )
      ),

      tabItem(
        tabName = "evolution",
        fluidRow(column(12, div(class = "section-header",
          h3("Network Evolution Over Time"),
          p("Tracking ln GDP per capita, openness, and network connectivity across European OECD countries.")
        ))),
        fluidRow(
          column(8, div(class = "viz-container", plotlyOutput("outcomeEvolutionPlot", height = "360px"))),
          column(4, div(class = "viz-container", h4("Connectivity Over Time"), plotlyOutput("connectivityPlot", height = "310px")))
        ),
        fluidRow(style = "margin-top: 15px;",
          column(12, div(class = "viz-container", h4("Selected Country GDP Trajectories"), plotlyOutput("countryTrajectoryPlot", height = "360px")))
        )
      ),

      tabItem(
        tabName = "analysis",
        fluidRow(column(12, div(class = "section-header",
          h3("CCE Control-Function and CFQR-CCE Analysis"),
          p("The auxiliary equation for openness is estimated by CCE-OLS. The outcome equation is estimated by unit-specific CFQR-CCE at tau = 0.5.")
        ))),
        fluidRow(
          column(8, div(class = "viz-container", plotlyOutput("opennessGdpScatter", height = "450px"))),
          column(4, div(class = "viz-container", h4("Top Countries by ln GDP per capita"), plotlyOutput("topCountriesPlot", height = "420px")))
        ),
        fluidRow(style = "margin-top: 15px;",
          column(12, div(class = "viz-container", h4("Estimated Main CFQR-CCE Coefficients"), plotlyOutput("coefPlot", height = "380px")))
        ),
        fluidRow(style = "margin-top: 15px;",
          column(12, div(class = "info-box-custom",
            h4("Model construction"),
            p(tags$b("Network support: "), "countries are potential neighbours when their longitude/latitude lattice distance satisfies d_ij <= d_c."),
            p(tags$b("Time-varying edge strength: "), "among geographically/topologically close countries, weights are stronger when observed openness levels are similar."),
            p(tags$b("Control function: "), "the residual epshat_it comes from a country-specific CCE-OLS regression of ln_openness on X_z and nodal averages."),
            p(tags$b("Outcome equation: "), "ln GDP per capita is estimated using CFQR-CCE with network terms, lagged GDP, X_y variables, CCE nodal averages, and epshat_it.")
          ))
        )
      )
    )
  )
)

############################################################
## Server
############################################################

server <- function(input, output, session) {

  current_index <- reactive({ match(input$year, years) })

  current_panel <- reactive({
    panel_plot %>% filter(year_actual == input$year)
  })

  current_summary <- reactive({
    network_summary %>% filter(year_actual == input$year)
  })

  output$avgYText <- renderText({ round(current_summary()$avg_Y, 3) })
  output$edgeCountText <- renderText({ current_summary()$edges })
  output$modelSizeText <- renderText({ paste0(N, " countries; pT = ", pT) })

  output$networkPlot <- renderPlotly({
    tt <- current_index()
    df_t <- current_panel() %>% arrange(match(country, countries))
    make_network_plotly(
      W = Wlist[[tt]],
      coords = coords,
      Y_t = as.numeric(Y_original_mat[, tt]),
      Z_t = as.numeric(Z_original_mat[, tt]),
      epshat_t = df_t$epshat,
      year_t = input$year,
      country_names = countries,
      dc = DC_DEGREES
    )
  })

  output$statCountries <- renderText({ N })
  output$statEdges <- renderText({ current_summary()$edges })
  output$statDensity <- renderText({ paste0(round(100 * current_summary()$density, 1), "%") })
  output$statAvgOutcome <- renderText({ round(current_summary()$avg_Y, 3) })

  output$outcomeEvolutionPlot <- renderPlotly({
    plot_ly(network_summary) %>%
      add_trace(
        x = ~year_actual, y = ~avg_Y, type = "scatter", mode = "lines+markers",
        name = "Average ln GDP per capita",
        line = list(width = 3, color = "rgba(44, 123, 182, 0.85)"),
        marker = list(size = 7, color = "rgba(44, 123, 182, 0.85)"),
        hovertemplate = "Year = %{x}<br>Average ln GDP per capita = %{y:.3f}<extra></extra>"
      ) %>%
      add_trace(
        x = ~year_actual, y = ~avg_Z, type = "scatter", mode = "lines+markers",
        name = "Average ln openness",
        line = list(width = 3, color = "rgba(215, 48, 39, 0.75)"),
        marker = list(size = 7, color = "rgba(215, 48, 39, 0.75)"),
        yaxis = "y2",
        hovertemplate = "Year = %{x}<br>Average ln openness = %{y:.3f}<extra></extra>"
      ) %>%
      layout(
        xaxis = list(title = "Year", dtick = 2, gridcolor = "rgba(230,230,230,0.7)"),
        yaxis = list(title = "Average ln GDP per capita", gridcolor = "rgba(230,230,230,0.7)"),
        yaxis2 = list(title = "Average ln openness", overlaying = "y", side = "right", showgrid = FALSE),
        legend = list(orientation = "h", x = 0.05, y = 1.12),
        margin = list(l = 80, r = 80, t = 30, b = 50),
        plot_bgcolor = "rgba(255,255,255,0.95)", paper_bgcolor = "white",
        font = list(family = "Gill Sans, Arial, sans-serif")
      ) %>% config(displayModeBar = FALSE)
  })

  output$connectivityPlot <- renderPlotly({
    plot_ly(network_summary) %>%
      add_trace(
        x = ~year_actual, y = ~edges, type = "bar",
        marker = list(color = "rgba(60, 141, 188, 0.75)", line = list(color = "rgba(255,255,255,0.8)", width = 1)),
        hovertemplate = "Year = %{x}<br>Edges = %{y}<extra></extra>"
      ) %>%
      layout(
        xaxis = list(title = "Year", dtick = 2, gridcolor = "rgba(230,230,230,0.7)"),
        yaxis = list(title = "Nonzero links", gridcolor = "rgba(230,230,230,0.7)"),
        margin = list(l = 60, r = 20, t = 20, b = 50),
        plot_bgcolor = "rgba(255,255,255,0.95)", paper_bgcolor = "white",
        font = list(family = "Gill Sans, Arial, sans-serif")
      ) %>% config(displayModeBar = FALSE)
  })

  output$countryTrajectoryPlot <- renderPlotly({
    selected <- intersect(c("Germany", "France", "United Kingdom", "Italy", "Spain", "Switzerland", "Sweden", "Norway"), countries)
    if (length(selected) < 4) selected <- head(countries, 6)

    df <- panel_plot %>% filter(country %in% selected)
    p <- plot_ly()
    for (cc in selected) {
      df_cc <- df %>% filter(country == cc)
      p <- p %>% add_trace(
        data = df_cc, x = ~year_actual, y = ~Y, type = "scatter", mode = "lines+markers",
        name = cc, line = list(width = 2), marker = list(size = 6),
        hovertemplate = paste0(cc, "<br>Year = %{x}<br>ln GDP per capita = %{y:.3f}<extra></extra>")
      )
    }
    p %>% layout(
      xaxis = list(title = "Year", dtick = 2, gridcolor = "rgba(230,230,230,0.7)"),
      yaxis = list(title = "ln GDP per capita", gridcolor = "rgba(230,230,230,0.7)"),
      legend = list(orientation = "h", x = 0.05, y = 1.12),
      margin = list(l = 70, r = 30, t = 30, b = 50),
      plot_bgcolor = "rgba(255,255,255,0.95)", paper_bgcolor = "white",
      font = list(family = "Gill Sans, Arial, sans-serif")
    ) %>% config(displayModeBar = FALSE)
  })

  output$opennessGdpScatter <- renderPlotly({
    df <- current_panel() %>% arrange(match(country, countries))
    W <- Wlist[[current_index()]]
    peer_gdp <- as.numeric(W %*% Y_original_mat[, current_index()])
    df <- df %>% mutate(peer_gdp = peer_gdp)

    plot_ly(df) %>%
      add_trace(
        x = ~Z, y = ~Y, type = "scatter", mode = "markers",
        text = ~paste0(
          country,
          "<br>ln openness Z = ", round(Z, 3),
          "<br>ln GDP per capita Y = ", round(Y, 3),
          "<br>Network GDP exposure WY = ", round(peer_gdp, 3),
          "<br>epshat from CCE-Z = ", round(epshat, 3)
        ),
        hoverinfo = "text",
        marker = list(
          size = ~scale_range(peer_gdp, 12, 42),
          color = ~peer_gdp,
          colorscale = list(c(0, "blue"), c(0.5, "yellow"), c(1, "red")),
          showscale = TRUE,
          colorbar = list(title = list(text = "Network GDP<br>exposure"), thickness = 14, len = 0.75, x = 1.02),
          line = list(color = "black", width = 1), opacity = 0.85
        )
      ) %>%
      layout(
        xaxis = list(title = "Auxiliary variable Z_it = ln openness", gridcolor = "rgba(230,230,230,0.7)"),
        yaxis = list(title = "Outcome Y_it = ln GDP per capita", gridcolor = "rgba(230,230,230,0.7)"),
        margin = list(l = 75, r = 110, t = 30, b = 70),
        plot_bgcolor = "rgba(255,255,255,0.95)", paper_bgcolor = "white",
        font = list(family = "Gill Sans, Arial, sans-serif")
      ) %>% config(displayModeBar = FALSE)
  })

  output$topCountriesPlot <- renderPlotly({
    df <- current_panel() %>% arrange(desc(Y)) %>% slice(1:10) %>% mutate(country = factor(country, levels = rev(country)))
    plot_ly(df) %>%
      add_trace(
        x = ~Y, y = ~country, type = "bar", orientation = "h",
        marker = list(color = "rgba(243, 156, 18, 0.75)", line = list(color = "rgba(255,255,255,0.8)", width = 1)),
        text = ~round(Y, 3), textposition = "outside",
        hovertemplate = "%{y}<br>ln GDP per capita = %{x:.3f}<extra></extra>"
      ) %>%
      layout(
        xaxis = list(title = "ln GDP per capita", gridcolor = "rgba(230,230,230,0.7)"),
        yaxis = list(title = ""),
        margin = list(l = 120, r = 50, t = 20, b = 50),
        plot_bgcolor = "rgba(255,255,255,0.95)", paper_bgcolor = "white",
        font = list(family = "Gill Sans, Arial, sans-serif")
      ) %>% config(displayModeBar = FALSE)
  })

  output$coefPlot <- renderPlotly({
    df <- coef_summary %>%
      pivot_longer(cols = c(Gamma_WY, lambda_Ylag, rho_WlagYlag, delta_epshat), names_to = "coefficient", values_to = "estimate") %>%
      filter(!is.na(estimate))

    plot_ly(df) %>%
      add_trace(
        x = ~estimate, y = ~coefficient, type = "box", orientation = "h",
        boxpoints = "all", jitter = 0.35, pointpos = 0,
        marker = list(color = "rgba(60,141,188,0.65)"),
        line = list(color = "rgba(44,62,80,0.85)"),
        hovertemplate = "Coefficient: %{y}<br>Estimate: %{x:.3f}<extra></extra>"
      ) %>%
      layout(
        xaxis = list(title = "Country-specific CFQR-CCE estimates", gridcolor = "rgba(230,230,230,0.7)"),
        yaxis = list(title = ""),
        margin = list(l = 140, r = 30, t = 30, b = 60),
        plot_bgcolor = "rgba(255,255,255,0.95)", paper_bgcolor = "white",
        font = list(family = "Gill Sans, Arial, sans-serif")
      ) %>% config(displayModeBar = FALSE)
  })
}

############################################################
## Run application
############################################################

shinyApp(ui, server)
