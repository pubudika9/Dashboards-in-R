############################################################
## University Student Social Network Dashboard
##
## Y_it    = education outcome of student i at time t
## Z_it    = engagement / social-academic activity index
## W_ij,t  = peer influence weight from student j to i
##
## Dashboard period: t = 0, ..., 15
############################################################

rm(list = ls())

############################################################
## Required packages
############################################################

library(shiny)
library(shinydashboard)
library(plotly)
library(dplyr)
library(tidyr)
library(igraph)
library(scales)

############################################################
## 1. Expanding-domain social-opportunity coordinates
############################################################

make_student_social_coordinates <- function(
    N,
    d0 = 1.25,
    jitter = 0.10,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  side <- ceiling(sqrt(N))
  
  grid <- expand.grid(
    academic_group = 0:(side - 1),
    social_group   = 0:(side - 1)
  )
  
  coords <- as.matrix(grid[1:N, ])
  
  coords <- d0 * coords
  
  coords <- coords + matrix(
    runif(2 * N, -jitter, jitter),
    nrow = N,
    ncol = 2
  )
  
  colnames(coords) <- c(
    "academic_proximity",
    "social_proximity"
  )
  
  return(coords)
}

############################################################
## 2. L-infinity distance matrix
############################################################

make_linf_distance_matrix <- function(coords) {
  N <- nrow(coords)
  
  D <- matrix(0, N, N)
  
  for (i in 1:N) {
    for (j in 1:N) {
      D[i, j] <- max(abs(coords[i, ] - coords[j, ]))
    }
  }
  
  return(D)
}

############################################################
## 3. Kernel h(Z_it, Z_jt)
############################################################

h_kernel <- function(z_i, z_j, eta = 1.0) {
  exp(-eta * abs(z_i - z_j))
}

############################################################
## 4. Build W_t according to Assumption 3
############################################################

build_Wt_assumption3 <- function(
    z_t,
    Dmat,
    dc = 2.75,
    eta = 1.0,
    weight_threshold = 0.02,
    renormalize_after_threshold = TRUE
) {
  N <- length(z_t)
  
  H <- outer(
    z_t,
    z_t,
    Vectorize(function(a, b) h_kernel(a, b, eta = eta))
  )
  
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
  
  return(W)
}

############################################################
## 5. Build W_t for all t
############################################################

build_Wlist_assumption3 <- function(
    Z,
    Dmat,
    dc = 2.75,
    eta = 1.0,
    weight_threshold = 0.02,
    renormalize_after_threshold = TRUE
) {
  Tt <- ncol(Z)
  
  Wlist <- vector("list", Tt)
  
  for (tt in 1:Tt) {
    Wlist[[tt]] <- build_Wt_assumption3(
      z_t = Z[, tt],
      Dmat = Dmat,
      dc = dc,
      eta = eta,
      weight_threshold = weight_threshold,
      renormalize_after_threshold = renormalize_after_threshold
    )
  }
  
  return(Wlist)
}

############################################################
## 6. Student education outcome DGP
############################################################

generate_student_social_network_data <- function(
    N,
    Tt,
    tau = 0.5,
    
    lambda_true = 0.25,
    gamma_true  = 0.35,
    rho_true    = 0.15,
    delta_true  = 0.60,
    
    phi_f = 0.70,
    
    d0 = 1.25,
    jitter = 0.10,
    dc = 2.75,
    
    eta = 1.0,
    weight_threshold = 0.02,
    seed_coords = 123
) {
  
  Sburn <- 200
  Tfull <- Tt + Sburn
  
  m <- 3
  py <- 2
  pz <- 3
  
  rho_x <- 0.8
  rho_eps <- 0.3
  rho_E <- 0.3
  coef_spread <- 0.04
  
  coords <- make_student_social_coordinates(
    N = N,
    d0 = d0,
    jitter = jitter,
    seed = seed_coords
  )
  
  Dmat <- make_linf_distance_matrix(coords)
  
  ##########################################################
  ## University-level latent factors
  ##########################################################
  
  Fmat <- matrix(0, nrow = Tfull, ncol = m)
  
  for (k in 1:m) {
    for (tt in 2:Tfull) {
      Fmat[tt, k] <- phi_f * Fmat[tt - 1, k] +
        sqrt(1 - phi_f^2) * rnorm(1)
    }
  }
  
  colnames(Fmat) <- c(
    "assessment_difficulty",
    "teaching_resource_shock",
    "campus_motivation_shock"
  )
  
  ##########################################################
  ## Student-specific loadings
  ##########################################################
  
  alpha_x <- matrix(0.5 + rnorm(N * py), nrow = N, ncol = py)
  
  gamma_x <- array(
    0.5 + rnorm(N * py * m),
    dim = c(N, py, m)
  )
  
  gamma_y <- matrix(0.5 + rnorm(N * m), nrow = N, ncol = m)
  gamma_z <- matrix(0.5 + rnorm(N * m), nrow = N, ncol = m)
  
  ##########################################################
  ## Outcome regressors
  ## Xy[,,1] = study effort
  ## Xy[,,2] = attendance / participation
  ##########################################################
  
  v <- array(0, dim = c(N, Tfull, py))
  
  for (k in 1:py) {
    for (tt in 2:Tfull) {
      v[, tt, k] <- rho_x * v[, tt - 1, k] +
        sqrt(1 - rho_x^2) * rnorm(N)
    }
  }
  
  Xy <- array(0, dim = c(N, Tfull, py))
  
  for (tt in 1:Tfull) {
    for (k in 1:py) {
      Xy[, tt, k] <- alpha_x[, k] +
        gamma_x[, k, ] %*% Fmat[tt, ] +
        v[, tt, k]
    }
  }
  
  study_effort <- Xy[, , 1]
  attendance_participation <- Xy[, , 2]
  
  ##########################################################
  ## Auxiliary regressors for Z_it
  ##########################################################
  
  campus_participation <- matrix(0, nrow = N, ncol = Tfull)
  online_learning_activity <- matrix(0, nrow = N, ncol = Tfull)
  
  for (tt in 2:Tfull) {
    campus_participation[, tt] <-
      0.4 * campus_participation[, tt - 1] +
      sqrt(1 - 0.4^2) * rnorm(N)
    
    online_learning_activity[, tt] <-
      0.4 * online_learning_activity[, tt - 1] +
      sqrt(1 - 0.4^2) * rnorm(N)
  }
  
  Xz_no_intercept <- array(0, dim = c(N, Tfull, pz))
  
  Xz_no_intercept[, , 1] <- campus_participation
  Xz_no_intercept[, , 2] <- online_learning_activity
  Xz_no_intercept[, , 3] <- study_effort
  
  beta_z_centre <- c(
    0.8,
    0.6,
    0.5
  )
  
  beta_z_i <- matrix(
    rep(beta_z_centre, each = N),
    nrow = N,
    ncol = pz
  ) + 0.10 * matrix(rnorm(N * pz), nrow = N, ncol = pz)
  
  ##########################################################
  ## Unobserved networking tendency / motivation
  ##########################################################
  
  eps_full <- matrix(0, nrow = N, ncol = Tfull)
  
  for (tt in 2:Tfull) {
    eps_full[, tt] <- rho_eps * eps_full[, tt - 1] +
      sqrt(1 - rho_eps^2) * rnorm(N)
  }
  
  ##########################################################
  ## Engagement / social-academic activity index Z_it
  ##########################################################
  
  Zfull <- matrix(0, nrow = N, ncol = Tfull)
  
  for (tt in 1:Tfull) {
    Xz_t <- Xz_no_intercept[, tt, ]
    
    Zfull[, tt] <- rowSums(beta_z_i * Xz_t) +
      gamma_z %*% Fmat[tt, ] +
      eps_full[, tt]
  }
  
  Wlist_full <- build_Wlist_assumption3(
    Z = Zfull,
    Dmat = Dmat,
    dc = dc,
    eta = eta,
    weight_threshold = weight_threshold,
    renormalize_after_threshold = TRUE
  )
  
  ##########################################################
  ## Structural coefficients
  ##########################################################
  
  Gamma_i_true <- lambda_true +
    coef_spread * (2 * runif(N) - 1)
  
  lambda_i_true <- gamma_true +
    coef_spread * (2 * runif(N) - 1)
  
  rho_i_true <- rho_true +
    coef_spread * (2 * runif(N) - 1)
  
  delta_i_true <- delta_true +
    coef_spread * (2 * runif(N) - 1)
  
  beta_study_effort <- 0.50 + 0.20 * (1:N) / N
  beta_attendance <- 0.40 + 0.20 * (1:N) / N
  
  keep_idx_full <- (Sburn + 1):Tfull
  
  fbar <- colMeans(Fmat[keep_idx_full, , drop = FALSE])
  effort_bar_i <- rowMeans(study_effort[, keep_idx_full, drop = FALSE])
  
  student_ability <- rnorm(N)
  
  alpha_i_true <- effort_bar_i +
    gamma_y %*% fbar +
    student_ability
  
  ##########################################################
  ## Quantile disturbance
  ##########################################################
  
  Eraw <- matrix(0, nrow = N, ncol = Tfull)
  
  for (tt in 2:Tfull) {
    Eraw[, tt] <- rho_E * Eraw[, tt - 1] +
      sqrt(1 - rho_E^2) * rnorm(N)
  }
  
  E_full <- Eraw - qnorm(tau)
  
  ##########################################################
  ## Generate education outcomes
  ##########################################################
  
  Yfull <- matrix(0, nrow = N, ncol = Tfull)
  
  for (tt in 1:Tfull) {
    
    if (tt == 1) {
      Y_prev <- rep(0, N)
      WlagY_prev <- rep(0, N)
    } else {
      Y_prev <- Yfull[, tt - 1]
      WlagY_prev <- Wlist_full[[tt - 1]] %*% Yfull[, tt - 1]
    }
    
    effort_t <- study_effort[, tt]
    attendance_t <- attendance_participation[, tt]
    
    common_t <- gamma_y %*% Fmat[tt, ]
    
    rhs_t <- alpha_i_true +
      lambda_i_true * Y_prev +
      rho_i_true * as.numeric(WlagY_prev) +
      beta_study_effort * effort_t +
      beta_attendance * attendance_t +
      as.numeric(common_t) +
      delta_i_true * eps_full[, tt] +
      E_full[, tt]
    
    A_t <- diag(N) - diag(Gamma_i_true) %*% Wlist_full[[tt]]
    
    Yfull[, tt] <- solve(A_t, rhs_t)
  }
  
  keep <- (Sburn + 1):Tfull
  
  Y <- Yfull[, keep, drop = FALSE]
  Z <- Zfull[, keep, drop = FALSE]
  Wlist <- Wlist_full[keep]
  
  X <- array(0, dim = c(N, Tt, 4))
  
  X[, , 1] <- 1
  X[, , 3] <- study_effort[, keep]
  X[, , 4] <- attendance_participation[, keep]
  
  X[, 1, 2] <- 0
  
  for (tt in 2:Tt) {
    X[, tt, 2] <- Y[, tt - 1]
  }
  
  Xz <- array(0, dim = c(N, Tt, pz + 1))
  Xz[, , 1] <- 1
  
  for (k in 1:pz) {
    Xz[, , k + 1] <- Xz_no_intercept[, keep, k]
  }
  
  return(
    list(
      Y = Y,
      Z = Z,
      X = X,
      Xz = Xz,
      Wlist = Wlist,
      coords = coords,
      Dmat = Dmat,
      F = Fmat[keep, , drop = FALSE],
      eps_true = eps_full[, keep, drop = FALSE],
      
      study_effort = study_effort[, keep, drop = FALSE],
      attendance_participation = attendance_participation[, keep, drop = FALSE],
      campus_participation = campus_participation[, keep, drop = FALSE],
      online_learning_activity = online_learning_activity[, keep, drop = FALSE],
      
      Gamma_i_true = Gamma_i_true,
      lambda_i_true = lambda_i_true,
      rho_i_true = rho_i_true,
      delta_i_true = delta_i_true,
      beta_study_effort = beta_study_effort,
      beta_attendance = beta_attendance,
      alpha_i_true = alpha_i_true,
      student_ability = student_ability
    )
  )
}

############################################################
## 7. Simulate dummy data for t = 0, ..., 15
############################################################

set.seed(123)

N <- 30
T_periods <- 16

dat <- generate_student_social_network_data(
  N = N,
  Tt = T_periods,
  tau = 0.5,
  lambda_true = 0.25,
  gamma_true = 0.35,
  rho_true = 0.15,
  delta_true = 0.60,
  phi_f = 0.70,
  d0 = 1.25,
  jitter = 0.10,
  dc = 2.75,
  eta = 1.0,
  weight_threshold = 0.02,
  seed_coords = 123
)

time_values <- 0:15

students <- data.frame(
  student_id = paste0("S", 1:N),
  i = 1:N,
  academic_proximity = dat$coords[, 1],
  social_proximity = dat$coords[, 2],
  ability = as.numeric(dat$student_ability),
  beta_study_effort = as.numeric(dat$beta_study_effort),
  beta_attendance = as.numeric(dat$beta_attendance)
)

panel_data <- expand.grid(
  i = 1:N,
  time = time_values
) %>%
  arrange(time, i) %>%
  mutate(
    student_id = paste0("S", i),
    Y = as.vector(dat$Y[, 1:T_periods]),
    Z = as.vector(dat$Z[, 1:T_periods]),
    study_effort = as.vector(dat$study_effort[, 1:T_periods]),
    attendance = as.vector(dat$attendance_participation[, 1:T_periods]),
    campus_participation = as.vector(dat$campus_participation[, 1:T_periods]),
    online_learning_activity = as.vector(dat$online_learning_activity[, 1:T_periods])
  ) %>%
  left_join(students, by = c("i", "student_id"))

network_summary <- data.frame(
  time = time_values,
  edges = sapply(dat$Wlist[1:T_periods], function(W) sum(W > 0)),
  density = sapply(dat$Wlist[1:T_periods], function(W) {
    sum(W > 0) / (N * (N - 1))
  }),
  avg_weight = sapply(dat$Wlist[1:T_periods], function(W) {
    ifelse(sum(W > 0) == 0, NA_real_, mean(W[W > 0]))
  }),
  max_weight = sapply(dat$Wlist[1:T_periods], function(W) max(W)),
  avg_Y = sapply(1:T_periods, function(tt) mean(dat$Y[, tt])),
  sd_Y = sapply(1:T_periods, function(tt) sd(dat$Y[, tt])),
  avg_Z = sapply(1:T_periods, function(tt) mean(dat$Z[, tt]))
)

############################################################
## 8. Helper functions for Plotly network visualization
############################################################

scale_range <- function(x, a, b) {
  x <- as.numeric(x)
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  
  if (abs(xmax - xmin) < 1e-12) {
    return(rep((a + b) / 2, length(x)))
  }
  
  a + (x - xmin) / (xmax - xmin) * (b - a)
}

make_network_plotly <- function(
    W,
    coords,
    Y_t,
    Z_t,
    time_t,
    show_labels = TRUE
) {
  
  N <- nrow(W)
  
  node_df <- data.frame(
    i = 1:N,
    student_id = paste0("S", 1:N),
    x = coords[, 1],
    y = coords[, 2],
    Y = Y_t,
    Z = Z_t,
    node_size = scale_range(Y_t, 14, 42)
  )
  
  edge_idx <- which(W > 0, arr.ind = TRUE)
  
  p <- plot_ly()
  
  if (nrow(edge_idx) > 0) {
    
    edge_weights <- W[edge_idx]
    edge_widths <- scale_range(edge_weights, 0.5, 5)
    
    for (k in seq_len(nrow(edge_idx))) {
      
      i <- edge_idx[k, 1]
      j <- edge_idx[k, 2]
      
      p <- p %>%
        add_trace(
          x = c(coords[i, 1], coords[j, 1], NA),
          y = c(coords[i, 2], coords[j, 2], NA),
          type = "scatter",
          mode = "lines",
          line = list(
            width = edge_widths[k],
            color = "rgba(80,80,80,0.30)"
          ),
          hoverinfo = "text",
          text = paste0(
            "Directed peer link: S", j, " → S", i,
            "<br>Weight W<sub>", i, j, ",", time_t, "</sub> = ",
            round(W[i, j], 3)
          ),
          showlegend = FALSE
        )
    }
  }
  
  p <- p %>%
    add_trace(
      data = node_df,
      x = ~x,
      y = ~y,
      type = "scatter",
      mode = ifelse(show_labels, "markers+text", "markers"),
      text = ~student_id,
      textposition = "middle center",
      textfont = list(size = 10, color = "black"),
      marker = list(
        size = ~node_size,
        color = ~Y,
        colorscale = list(
          c(0, "blue"),
          c(0.2, "cyan"),
          c(0.4, "green"),
          c(0.6, "yellow"),
          c(0.8, "orange"),
          c(1, "red")
        ),
        showscale = TRUE,
        colorbar = list(
          title = list(text = "Education<br>outcome"),
          thickness = 14,
          len = 0.75,
          x = 1.02
        ),
        line = list(color = "black", width = 1)
      ),
      hovertemplate = paste0(
        "<b>%{text}</b>",
        "<br>Education outcome Y: %{marker.color:.2f}",
        "<br>Academic proximity: %{x:.2f}",
        "<br>Social proximity: %{y:.2f}",
        "<extra></extra>"
      ),
      showlegend = FALSE
    ) %>%
    layout(
      title = list(
        text = paste0("University Student Social Network at t = ", time_t),
        x = 0.02
      ),
      xaxis = list(
        title = "Academic proximity coordinate",
        zeroline = FALSE,
        gridcolor = "rgba(220,220,220,0.6)"
      ),
      yaxis = list(
        title = "Social proximity coordinate",
        zeroline = FALSE,
        gridcolor = "rgba(220,220,220,0.6)",
        scaleanchor = "x",
        scaleratio = 1
      ),
      plot_bgcolor = "rgba(255,255,255,0.95)",
      paper_bgcolor = "white",
      margin = list(l = 60, r = 110, t = 60, b = 60),
      font = list(family = "Gill Sans, Arial, sans-serif")
    ) %>%
    config(displayModeBar = FALSE)
  
  return(p)
}

############################################################
## 9. Shiny UI
############################################################

ui <- dashboardPage(
  title = "Student Social Network",
  skin = "black",
  
  dashboardHeader(
    title = "Student Social Network Dashboard",
    titleWidth = 280
  ),
  
  dashboardSidebar(
    width = 280,
    
    sidebarMenu(
      id = "tabs",
      menuItem("Current Network", tabName = "current", icon = icon("project-diagram")),
      menuItem("Evolution", tabName = "evolution", icon = icon("chart-line")),
      menuItem("Peer Influence Analysis", tabName = "impact", icon = icon("users"))
    ),
    
    div(
      style = "padding: 15px; color: #d8d8d8;",
      p("Explore how a simulated university student social network evolves over time."),
      p("Node size and colour represent education outcome. Edge width represents peer influence weight.")
    ),
    
    div(
      style = "padding: 15px;",
      h5("Select Time Period"),
      sliderInput(
        "time",
        NULL,
        min = 0,
        max = 15,
        value = 0,
        step = 1,
        sep = "",
        animate = FALSE,
        width = "100%",
        ticks = TRUE
      ),
      actionButton(
        "playBtn",
        label = "",
        icon = icon("play"),
        class = "btn btn-primary"
      )
    ),
    
    div(
      style = "padding: 15px;",
      h5("Network Statistics"),
      strong("Average education outcome"),
      div(style = "font-size: 26px; color: #00c0ef;", textOutput("avgYText")),
      br(),
      strong("Nonzero peer links"),
      div(style = "font-size: 26px; color: #f39c12;", textOutput("edgeCountText"))
    ),
    
    div(
      style = "position:absolute; bottom:20px; padding: 15px; color:#b8b8b8;",
      p("Data: Simulated student network"),
      p("Model: Endogenous network with Assumption 3")
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        * {
          font-family: 'Gill Sans', Arial, sans-serif !important;
        }

        body, .content-wrapper {
          background-color: #f7f7f7;
        }

        .section-header {
          background: white;
          padding: 18px 22px;
          margin-bottom: 15px;
          border-radius: 8px;
          box-shadow: 0 1px 5px rgba(0,0,0,0.08);
        }

        .section-subtitle {
          color: #666;
          font-size: 15px;
          margin-top: 5px;
        }

        .viz-container {
          background: white;
          padding: 15px;
          border-radius: 8px;
          box-shadow: 0 1px 5px rgba(0,0,0,0.08);
        }

        .stat-box {
          background: white;
          text-align: center;
          padding: 18px;
          border-radius: 8px;
          box-shadow: 0 1px 5px rgba(0,0,0,0.08);
        }

        .stat-value {
          font-size: 30px;
          font-weight: bold;
          color: #2c3e50;
        }

        .stat-label {
          font-size: 13px;
          color: #777;
          margin-top: 5px;
        }

        .info-box-custom {
          background: white;
          padding: 15px;
          border-left: 4px solid #3c8dbc;
          border-radius: 6px;
          box-shadow: 0 1px 5px rgba(0,0,0,0.08);
          margin-bottom: 15px;
        }
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

              var slider = $('#time').data('ionRangeSlider');
              if (!slider) return;

              var currentVal = slider.result.from;
              var maxVal = slider.result.max;
              var minVal = slider.result.min;

              if (currentVal >= maxVal) {
                slider.update({from: minVal});
              }

              interval = setInterval(function() {
                var current = slider.result.from;

                if (current < maxVal) {
                  slider.update({from: current + 1});
                } else {
                  clearInterval(interval);
                  playing = false;
                  $('#playBtn').find('i').removeClass('fa-pause').addClass('fa-play');
                }
              }, 1500);

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
      
      ######################################################
      ## Current Network Tab
      ######################################################
      
      tabItem(
        tabName = "current",
        
        fluidRow(
          column(
            12,
            div(
              class = "section-header",
              h3("Current Student Peer Network"),
              div(
                class = "section-subtitle",
                "Education outcomes and time-varying peer influence links among university students"
              ),
              p(
                "Each node is a student. Node colour and size represent the education outcome ",
                tags$b("Y_it"),
                ". Edge width represents the peer influence weight ",
                tags$b("W_ij,t"),
                ". Links are generated from engagement similarity and social-academic proximity."
              )
            )
          )
        ),
        
        fluidRow(
          column(
            12,
            div(
              class = "viz-container",
              plotlyOutput("networkPlot", height = "650px")
            )
          )
        ),
        
        fluidRow(
          style = "margin-top: 15px;",
          column(
            3,
            div(
              class = "stat-box",
              div(class = "stat-value", textOutput("statStudents")),
              div(class = "stat-label", "Students")
            )
          ),
          column(
            3,
            div(
              class = "stat-box",
              div(class = "stat-value", textOutput("statEdges")),
              div(class = "stat-label", "Nonzero Directed Peer Links")
            )
          ),
          column(
            3,
            div(
              class = "stat-box",
              div(class = "stat-value", textOutput("statDensity")),
              div(class = "stat-label", "Network Density")
            )
          ),
          column(
            3,
            div(
              class = "stat-box",
              div(class = "stat-value", textOutput("statAvgOutcome")),
              div(class = "stat-label", "Average Education Outcome")
            )
          )
        )
      ),
      
      ######################################################
      ## Evolution Tab
      ######################################################
      
      tabItem(
        tabName = "evolution",
        
        fluidRow(
          column(
            12,
            div(
              class = "section-header",
              h3("Network Evolution Over Time"),
              p("Tracking education outcomes and peer-network structure from t = 0 to t = 15."),
              p(
                "The time series plots summarize how average academic performance, engagement, and network connectivity evolve."
              )
            )
          )
        ),
        
        fluidRow(
          column(
            8,
            div(
              class = "viz-container",
              plotlyOutput("outcomeEvolutionPlot", height = "360px")
            )
          ),
          column(
            4,
            div(
              class = "viz-container",
              h4("Connectivity Over Time"),
              plotlyOutput("connectivityPlot", height = "310px")
            )
          )
        ),
        
        fluidRow(
          style = "margin-top: 15px;",
          column(
            12,
            div(
              class = "viz-container",
              h4("Student-Level Education Outcome Trajectories"),
              plotlyOutput("studentTrajectoryPlot", height = "360px")
            )
          )
        )
      ),
      
      ######################################################
      ## Impact / Peer Influence Tab
      ######################################################
      
      tabItem(
        tabName = "impact",
        
        fluidRow(
          column(
            12,
            div(
              class = "section-header",
              h3("Peer Influence and Education Outcomes"),
              p(
                "This section explores the relationship between a student's education outcome, engagement index, and peer exposure."
              )
            )
          )
        ),
        
        fluidRow(
          column(
            8,
            div(
              class = "viz-container",
              plotlyOutput("engagementOutcomeScatter", height = "450px")
            )
          ),
          column(
            4,
            div(
              class = "viz-container",
              h4("Top Students by Education Outcome"),
              plotlyOutput("topStudentsPlot", height = "420px")
            )
          )
        ),
        
        fluidRow(
          style = "margin-top: 15px;",
          column(
            12,
            div(
              class = "info-box-custom",
              h4("Interpretation"),
              p(
                tags$b("Network support: "),
                "students are potential peers only if they are close in the social-academic opportunity space."
              ),
              p(
                tags$b("Time-varying edge strength: "),
                "among nearby students, links are stronger when engagement indices are more similar."
              ),
              p(
                tags$b("Education outcome: "),
                "the simulated outcome depends on own lagged performance, contemporaneous peer outcomes, lagged peer outcomes, study effort, attendance, common university shocks, and unobserved engagement/motivation."
              )
            )
          )
        )
      )
    )
  )
)

############################################################
## 10. Server
############################################################

server <- function(input, output, session) {
  
  current_index <- reactive({
    input$time + 1
  })
  
  current_panel <- reactive({
    panel_data %>% filter(time == input$time)
  })
  
  current_W <- reactive({
    dat$Wlist[[current_index()]]
  })
  
  current_summary <- reactive({
    network_summary %>% filter(time == input$time)
  })
  
  ##########################################################
  ## Sidebar statistics
  ##########################################################
  
  output$avgYText <- renderText({
    round(current_summary()$avg_Y, 2)
  })
  
  output$edgeCountText <- renderText({
    current_summary()$edges
  })
  
  ##########################################################
  ## Current Network Tab
  ##########################################################
  
  output$networkPlot <- renderPlotly({
    tt <- current_index()
    
    make_network_plotly(
      W = dat$Wlist[[tt]],
      coords = dat$coords,
      Y_t = dat$Y[, tt],
      Z_t = dat$Z[, tt],
      time_t = input$time,
      show_labels = TRUE
    )
  })
  
  output$statStudents <- renderText({
    N
  })
  
  output$statEdges <- renderText({
    current_summary()$edges
  })
  
  output$statDensity <- renderText({
    paste0(round(100 * current_summary()$density, 1), "%")
  })
  
  output$statAvgOutcome <- renderText({
    round(current_summary()$avg_Y, 2)
  })
  
  ##########################################################
  ## Evolution Tab
  ##########################################################
  
  output$outcomeEvolutionPlot <- renderPlotly({
    
    plot_ly(network_summary) %>%
      add_trace(
        x = ~time,
        y = ~avg_Y,
        type = "scatter",
        mode = "lines+markers",
        name = "Average education outcome",
        line = list(width = 3, color = "rgba(44, 123, 182, 0.85)"),
        marker = list(size = 8, color = "rgba(44, 123, 182, 0.85)"),
        hovertemplate = "t = %{x}<br>Average Y = %{y:.2f}<extra></extra>"
      ) %>%
      add_trace(
        x = ~time,
        y = ~avg_Z,
        type = "scatter",
        mode = "lines+markers",
        name = "Average engagement index",
        line = list(width = 3, color = "rgba(215, 48, 39, 0.75)"),
        marker = list(size = 8, color = "rgba(215, 48, 39, 0.75)"),
        yaxis = "y2",
        hovertemplate = "t = %{x}<br>Average Z = %{y:.2f}<extra></extra>"
      ) %>%
      layout(
        xaxis = list(
          title = "Time period",
          dtick = 1,
          gridcolor = "rgba(230,230,230,0.7)"
        ),
        yaxis = list(
          title = "Average education outcome",
          gridcolor = "rgba(230,230,230,0.7)"
        ),
        yaxis2 = list(
          title = "Average engagement index",
          overlaying = "y",
          side = "right",
          showgrid = FALSE
        ),
        legend = list(
          orientation = "h",
          x = 0.05,
          y = 1.12
        ),
        margin = list(l = 70, r = 70, t = 30, b = 50),
        plot_bgcolor = "rgba(255,255,255,0.95)",
        paper_bgcolor = "white",
        font = list(family = "Gill Sans, Arial, sans-serif")
      ) %>%
      config(displayModeBar = FALSE)
  })
  
  output$connectivityPlot <- renderPlotly({
    
    plot_ly(network_summary) %>%
      add_trace(
        x = ~time,
        y = ~edges,
        type = "bar",
        marker = list(
          color = "rgba(60, 141, 188, 0.75)",
          line = list(color = "rgba(255,255,255,0.8)", width = 1)
        ),
        hovertemplate = "t = %{x}<br>Edges = %{y}<extra></extra>"
      ) %>%
      layout(
        xaxis = list(
          title = "Time period",
          dtick = 1,
          gridcolor = "rgba(230,230,230,0.7)"
        ),
        yaxis = list(
          title = "Nonzero peer links",
          gridcolor = "rgba(230,230,230,0.7)"
        ),
        margin = list(l = 60, r = 20, t = 20, b = 50),
        plot_bgcolor = "rgba(255,255,255,0.95)",
        paper_bgcolor = "white",
        font = list(family = "Gill Sans, Arial, sans-serif")
      ) %>%
      config(displayModeBar = FALSE)
  })
  
  output$studentTrajectoryPlot <- renderPlotly({
    
    selected_students <- paste0("S", c(1, 6, 12, 18, 24, 30))
    
    df <- panel_data %>%
      filter(student_id %in% selected_students)
    
    p <- plot_ly()
    
    for (sid in selected_students) {
      df_s <- df %>% filter(student_id == sid)
      
      p <- p %>%
        add_trace(
          data = df_s,
          x = ~time,
          y = ~Y,
          type = "scatter",
          mode = "lines+markers",
          name = sid,
          line = list(width = 2),
          marker = list(size = 6),
          hovertemplate = paste0(
            sid,
            "<br>t = %{x}",
            "<br>Y = %{y:.2f}",
            "<extra></extra>"
          )
        )
    }
    
    p %>%
      layout(
        xaxis = list(
          title = "Time period",
          dtick = 1,
          gridcolor = "rgba(230,230,230,0.7)"
        ),
        yaxis = list(
          title = "Education outcome",
          gridcolor = "rgba(230,230,230,0.7)"
        ),
        legend = list(
          orientation = "h",
          x = 0.05,
          y = 1.12
        ),
        margin = list(l = 70, r = 30, t = 30, b = 50),
        plot_bgcolor = "rgba(255,255,255,0.95)",
        paper_bgcolor = "white",
        font = list(family = "Gill Sans, Arial, sans-serif")
      ) %>%
      config(displayModeBar = FALSE)
  })
  
  ##########################################################
  ## Impact / Peer Influence Tab
  ##########################################################
  
  output$engagementOutcomeScatter <- renderPlotly({
    
    df <- current_panel()
    
    W <- current_W()
    
    peer_exposure <- as.numeric(W %*% dat$Y[, current_index()])
    
    df <- df %>%
      mutate(
        peer_exposure = peer_exposure
      )
    
    plot_ly(df) %>%
      add_trace(
        x = ~Z,
        y = ~Y,
        type = "scatter",
        mode = "markers",
        text = ~paste0(
          student_id,
          "<br>Engagement index Z = ", round(Z, 2),
          "<br>Education outcome Y = ", round(Y, 2),
          "<br>Peer exposure WY = ", round(peer_exposure, 2),
          "<br>Study effort = ", round(study_effort, 2),
          "<br>Attendance = ", round(attendance, 2)
        ),
        hoverinfo = "text",
        marker = list(
          size = ~scale_range(peer_exposure, 12, 42),
          color = ~peer_exposure,
          colorscale = list(
            c(0, "blue"),
            c(0.5, "yellow"),
            c(1, "red")
          ),
          showscale = TRUE,
          colorbar = list(
            title = list(text = "Peer<br>exposure"),
            thickness = 14,
            len = 0.75,
            x = 1.02
          ),
          line = list(color = "black", width = 1),
          opacity = 0.85
        )
      ) %>%
      layout(
        xaxis = list(
          title = "Engagement / social-academic activity index Z_it",
          gridcolor = "rgba(230,230,230,0.7)"
        ),
        yaxis = list(
          title = "Education outcome Y_it",
          gridcolor = "rgba(230,230,230,0.7)"
        ),
        margin = list(l = 75, r = 100, t = 30, b = 70),
        plot_bgcolor = "rgba(255,255,255,0.95)",
        paper_bgcolor = "white",
        font = list(family = "Gill Sans, Arial, sans-serif")
      ) %>%
      config(displayModeBar = FALSE)
  })
  
  output$topStudentsPlot <- renderPlotly({
    
    df <- current_panel() %>%
      arrange(desc(Y)) %>%
      slice(1:10) %>%
      mutate(student_id = factor(student_id, levels = rev(student_id)))
    
    plot_ly(df) %>%
      add_trace(
        x = ~Y,
        y = ~student_id,
        type = "bar",
        orientation = "h",
        marker = list(
          color = "rgba(243, 156, 18, 0.75)",
          line = list(color = "rgba(255,255,255,0.8)", width = 1)
        ),
        text = ~round(Y, 2),
        textposition = "outside",
        hovertemplate = "%{y}<br>Education outcome = %{x:.2f}<extra></extra>"
      ) %>%
      layout(
        xaxis = list(
          title = "Education outcome",
          gridcolor = "rgba(230,230,230,0.7)"
        ),
        yaxis = list(
          title = ""
        ),
        margin = list(l = 70, r = 50, t = 20, b = 50),
        plot_bgcolor = "rgba(255,255,255,0.95)",
        paper_bgcolor = "white",
        font = list(family = "Gill Sans, Arial, sans-serif")
      ) %>%
      config(displayModeBar = FALSE)
  })
}

############################################################
## 11. Run application
############################################################

shinyApp(ui, server)