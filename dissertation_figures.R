# ============================================================================
# dissertation_figures.R
# Reproduces all 7 figures in "Model and Calibration-Window Risk in UK
# Longevity" from the raw HMD GBR_NP files.
#
#   Fig 3.1  mortality surface ln m(x,t) + period e65          (data only)
#   Fig 3.2  rates by age, selected years                      (data only)
#   Fig 3.3  rates over time at ages 65/75/85                  (data only)
#   Fig 3.4  annual improvement, ages 60-89, with means        (data only)
#   Fig 5.1  out-of-sample validation at age 75 (LC & CBD)     (StMoMo)
#   Fig 5.2  best-estimate annuity by model x window           (StMoMo, ~min)
#   Fig 5.3  longevity capital SCR% heatmap                    (StMoMo, ~min)
#
# Usage: place this file in the same directory as GBR_data/ and run it
# section by section (RStudio: Ctrl+Enter). Part 0 must be run first;
# Part 1 needs only the data; Parts 2-3 need StMoMo (Part 3 takes ~2-5 min).
# Output: 7 PNG files in ./figures/ plus capital_grid_R.csv (use it to
# update Table 5.2).
# Note: the numbers produced by Part 3 are the results of YOUR run; the
# figures in Tables 5.1/5.2/5.3 and in the text of the dissertation should
# be updated to match them. Small differences from any earlier prototype
# results are normal and genuine.
# ============================================================================

## ---- PART 0 . packages, palette, data -------------------------------------
# install.packages(c("ggplot2", "patchwork", "StMoMo"))   # install once
library(ggplot2)
library(patchwork)

DATA_DIR <- "GBR_data/STATS"
dir.create("figures", showWarnings = FALSE)
AGES  <- 55:89
YEARS <- 1961:2018

NAVY <- "#1B2A4A"; TEAL <- "#2A9D8F"; AMBER <- "#E9A23B"; ROSE <- "#C9596B"
MUTED <- "#5A6B7B"; GRID <- "#D8E0E8"; GREY <- "#9AA7B2"

theme_diss <- function(base = 12) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = GRID, linewidth = 0.4),
          axis.text  = element_text(colour = MUTED),
          axis.title = element_text(colour = NAVY),
          plot.title = element_text(colour = NAVY, face = "bold", size = base + 1),
          legend.text  = element_text(colour = MUTED),
          legend.title = element_text(colour = NAVY))
}

# -- Parse an HMD 1x1 file: skip the 2 header lines; Age contains "110+" --
read_hmd <- function(fname) {
  df <- read.table(file.path(DATA_DIR, fname), skip = 2, header = TRUE,
                   na.strings = ".", stringsAsFactors = FALSE)
  df$Age <- as.integer(sub("\\+", "", df$Age))
  df
}
# Pivot to an ages x years matrix (males)
hmd_matrix <- function(df, col = "Male", ages = AGES, years = YEARS) {
  sub <- df[df$Age %in% ages & df$Year %in% years, ]
  m <- tapply(sub[[col]], list(sub$Age, sub$Year), function(z) z[1])
  m[as.character(ages), as.character(years), drop = FALSE]
}

Dxt <- hmd_matrix(read_hmd("Deaths_1x1.txt"))      # death counts
Ext <- hmd_matrix(read_hmd("Exposures_1x1.txt"))   # central exposures
mxt <- Dxt / Ext                                   # central death rates
lt  <- read_hmd("mltper_1x1.txt")                  # male period life table (for e65)
stopifnot(!anyNA(mxt))                             # the window should be complete

## ---- PART 1 . Chapter 3 figures (data only) --------------------------------

# ---- Figure 3.1: ln m(x,t) surface + period e65 ----
surf <- expand.grid(Age = AGES, Year = YEARS)
surf$logm <- as.vector(log(mxt))
p_surf <- ggplot(surf, aes(Year, Age, fill = logm)) +
  geom_raster() +
  scale_fill_viridis_c(name = expression(ln~m[x*","*t])) +
  labs(title = "Log central death-rate surface", x = "Calendar year  t", y = "Age  x") +
  coord_cartesian(expand = FALSE) + theme_diss()

e65 <- lt[lt$Age == 65 & lt$Year %in% YEARS, c("Year", "ex")]
p_e65 <- ggplot(e65, aes(Year, ex)) +
  geom_line(colour = NAVY, linewidth = 1.1) +
  geom_vline(xintercept = 2011, colour = MUTED, linetype = "dotted") +
  annotate("text", x = 2011.6, y = min(e65$ex) + 0.5, label = "2011",
           colour = MUTED, size = 3, hjust = 0) +
  labs(title = "Period life expectancy at 65",
       x = "Calendar year  t", y = expression(e[65])) + theme_diss()

fig31 <- p_surf + p_e65 + plot_layout(widths = c(1.5, 1))
ggsave("figures/fig3_1_surface_e65.png", fig31, width = 10, height = 3.6, dpi = 200)

# ---- Figure 3.2: rates by age, selected years ----
sel_years <- c(1961, 1980, 2000, 2018)
d32 <- expand.grid(Age = AGES, Year = sel_years)
d32$m <- mapply(function(a, y) mxt[as.character(a), as.character(y)], d32$Age, d32$Year)
d32$Year <- factor(d32$Year)
fig32 <- ggplot(d32, aes(Age, m, colour = Year)) +
  geom_line(linewidth = 1.2) +
  scale_y_log10() +
  scale_colour_manual(values = c(GREY, AMBER, TEAL, NAVY), name = "year") +
  labs(title = "UK male mortality by age, selected years",
       x = "Age  x", y = expression(central~death~rate~~m[x]~~(log~scale))) +
  theme_diss()
ggsave("figures/fig3_2_rates_by_age.png", fig32, width = 6.6, height = 4.3, dpi = 200)

# ---- Figure 3.3: rates over time at fixed ages ----
sel_ages <- c(65, 75, 85)
d33 <- expand.grid(Year = YEARS, Age = sel_ages)
d33$m <- mapply(function(a, y) mxt[as.character(a), as.character(y)], d33$Age, d33$Year)
d33$Age <- factor(d33$Age, labels = paste("age", sel_ages))
fig33 <- ggplot(d33, aes(Year, m, colour = Age)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 2011, colour = MUTED, linetype = "dotted") +
  scale_y_log10() +
  scale_colour_manual(values = c(TEAL, AMBER, ROSE), name = NULL) +
  labs(title = "UK male mortality over time, fixed ages",
       x = "Calendar year  t", y = expression(m[x*","*t]~~(log~scale))) +
  theme_diss()
ggsave("figures/fig3_3_rates_by_year.png", fig33, width = 6.6, height = 4.3, dpi = 200)

# ---- Figure 3.4: annual improvement, ages 60-89 ----
# Annual improvement (%) = -(change in ln m), averaged over ages 60-89,
# matching the definition used in the dissertation.
imp_ages <- as.character(60:89)
lm_mat   <- log(mxt[imp_ages, ])
impr <- -100 * colMeans(lm_mat[, -1] - lm_mat[, -ncol(lm_mat)])   # 1962..2018
yr2  <- YEARS[-1]
roll3 <- as.numeric(stats::filter(impr, rep(1/3, 3), sides = 2))  # 3-year average
m1 <- mean(impr[yr2 >= 1991 & yr2 <= 2011])    # should be ~2.8-2.9%
m2 <- mean(impr[yr2 >= 2012 & yr2 <= 2018])    # should be ~0.8%
d34 <- data.frame(Year = yr2, impr = impr, roll = roll3)
fig34 <- ggplot(d34, aes(Year, impr)) +
  geom_col(fill = GRID, width = 0.9) +
  geom_line(aes(y = roll), colour = NAVY, linewidth = 1.1, na.rm = TRUE) +
  annotate("segment", x = 1991, xend = 2011, y = m1, yend = m1, colour = TEAL, linewidth = 1.6) +
  annotate("segment", x = 2012, xend = 2018, y = m2, yend = m2, colour = ROSE, linewidth = 1.6) +
  annotate("text", x = 1996, y = m1 + 0.9, colour = TEAL, size = 3.2,
           label = sprintf("1991\u20132011 mean %.1f%%", m1)) +
  annotate("text", x = 2010, y = m2 - 1.1, colour = ROSE, size = 3.2,
           label = sprintf("2011\u20132018 mean %.1f%%", m2)) +
  geom_vline(xintercept = 2011, colour = MUTED, linetype = "dotted") +
  geom_hline(yintercept = 0, colour = MUTED, linewidth = 0.4) +
  labs(title = "UK male mortality improvement, ages 60\u201389",
       x = "Calendar year  t", y = "annual improvement in m (%)") +
  theme_diss()
ggsave("figures/fig3_4_improvement.png", fig34, width = 6.6, height = 4.3, dpi = 200)

## ---- PART 2 . Figure 5.1: out-of-sample validation at age 75 ---------------
library(StMoMo)
NSIM <- 1000

fit_window <- function(model, yrs, use_initial = FALSE) {
  ys  <- as.character(yrs)
  E   <- if (use_initial) Ext[, ys] + 0.5 * Dxt[, ys] else Ext[, ys]   # CBD uses initial exposures
  wxt <- genWeightMat(ages = AGES, years = yrs, clip = 3)              # clip sparse corner cohorts
  fit(model, Dxt = Dxt[, ys], Ext = E, ages = AGES, years = yrs, wxt = wxt)
}

tr_years <- 1961:2000; va_years <- 2001:2018; H <- length(va_years)
LCfit_v  <- fit_window(lc(link = "log"),    tr_years)
CBDfit_v <- fit_window(cbd(link = "logit"), tr_years, use_initial = TRUE)

set.seed(2025); simLC  <- simulate(LCfit_v,  nsim = NSIM, h = H)
set.seed(2025); simCBD <- simulate(CBDfit_v, nsim = NSIM, h = H)
fLC  <- forecast(LCfit_v,  h = H)
fCBD <- forecast(CBDfit_v, h = H)

a75 <- as.character(75)
logm_band <- function(simRates, central, is_q = FALSE) {
  r  <- simRates[a75, , ]                       # h x nsim
  cc <- central[a75, ]
  if (is_q) { r <- -log(1 - r); cc <- -log(1 - cc) }   # CBD: convert q to m
  data.frame(Year = va_years, mid = log(cc),
             lo = log(apply(r, 1, quantile, 0.025)),
             hi = log(apply(r, 1, quantile, 0.975)))
}
bLC  <- logm_band(simLC$rates,  fLC$rates)
bCBD <- logm_band(simCBD$rates, fCBD$rates, is_q = TRUE)
obs75 <- data.frame(Year = va_years, logm = log(mxt[a75, as.character(va_years)]))

fig51 <- ggplot() +
  geom_ribbon(data = bLC,  aes(Year, ymin = lo, ymax = hi), fill = AMBER, alpha = 0.15) +
  geom_ribbon(data = bCBD, aes(Year, ymin = lo, ymax = hi), fill = TEAL,  alpha = 0.15) +
  geom_line(data = bLC,  aes(Year, mid, colour = "LC forecast"),  linewidth = 1.1) +
  geom_line(data = bCBD, aes(Year, mid, colour = "CBD forecast"), linewidth = 1.1) +
  geom_point(data = obs75, aes(Year, logm, colour = "actual"), size = 2) +
  geom_vline(xintercept = 2011, colour = MUTED, linetype = "dotted") +
  scale_colour_manual(NULL, values = c("actual" = NAVY, "LC forecast" = AMBER,
                                       "CBD forecast" = TEAL)) +
  labs(title = "Out-of-sample validation at age 75\ntrained on 1961\u20132000, predicting 2001\u20132018",
       x = "Calendar year  t", y = expression(ln~m[75*","*t])) +
  theme_diss()
ggsave("figures/fig5_1_validation.png", fig51, width = 6.6, height = 4.2, dpi = 200)

## ---- PART 3 . Figures 5.2 & 5.3: model x window -> annuity & capital -------
# As specified in Chapter 4 of the dissertation: male aged 65 in 2018,
# discount rate i = 1%, Gompertz tail to close the life table,
# SCR = VaR99.5%(L) - BEL, reported as a percentage of BEL.
# Runtime roughly 2-5 minutes.

windows <- list(Long = 1961:2018, Recent = 1991:2018, PreSlow = 1980:2010)
x0 <- 65; base_year <- 2018; top_age <- 110; v <- 1 / 1.01

# Gompertz extrapolation: regress log m on age over ages 75-89 and
# extrapolate to age > 89. Closed-form least squares, roughly an order of
# magnitude faster than lm() -- this is called ~180k times in the
# simulation loop.
gompertz_m <- function(mcol, age) {
  xs <- AGES[AGES >= 75]; ys <- log(mcol[AGES >= 75])
  b2 <- sum((xs - mean(xs)) * (ys - mean(ys))) / sum((xs - mean(xs))^2)
  b1 <- mean(ys) - b2 * mean(xs)
  exp(b1 + b2 * age)
}

# Look up m(age, year) from observed rates (inside the window) combined
# with forecast rates (after the window closes).
make_lookup <- function(win, fut_years, rate_mat) {
  # rate_mat: ages x length(fut_years) matrix of (simulated or point)
  # forecast central death rates
  obs_years <- win
  function(age, year) {
    col <- if (year <= max(obs_years)) mxt[, as.character(min(year, max(obs_years)))]
           else rate_mat[, as.character(min(year, max(fut_years)))]
    if (age <= 89) col[as.character(age)] else gompertz_m(col, age)
  }
}

annuity_e65 <- function(mfun) {
  kp <- numeric(top_age - x0); kp_now <- 1; ann <- 1; e <- 0   # k = 0 term: v^0 * 1
  for (k in 1:(top_age - x0 - 1)) {
    m <- mfun(x0 + k - 1, base_year + k - 1)
    kp_now <- kp_now * exp(-m)                                  # survival p = exp(-m)
    ann <- ann + v^k * kp_now
    e   <- e + kp_now
  }
  c(annuity = ann, e65 = e + 0.5)
}

run_cell <- function(model_name, win) {
  endY <- max(win); Hc <- (base_year + (top_age - x0)) - endY
  fut  <- endY + seq_len(Hc)
  is_q <- model_name == "CBD"
  fitc <- switch(model_name,
    LC  = fit_window(lc(link = "log"), win),
    CBD = fit_window(cbd(link = "logit"), win, use_initial = TRUE),
    RH  = { LC0 <- fit_window(lc(link = "log"), win)
            ys <- as.character(win)
            fit(rh(link = "log", cohortAgeFun = "1", approxConst = TRUE),
                Dxt = Dxt[, ys], Ext = Ext[, ys], ages = AGES, years = win,
                wxt = genWeightMat(AGES, win, clip = 3),
                start.ax = LC0$ax, start.bx = LC0$bx, start.kt = LC0$kt) })
  fc <- if (model_name == "RH") forecast(fitc, h = Hc, gc.order = c(1, 1, 0))
        else                    forecast(fitc, h = Hc)
  set.seed(2025)
  sm <- if (model_name == "RH") simulate(fitc, nsim = NSIM, h = Hc, gc.order = c(1, 1, 0))
        else                    simulate(fitc, nsim = NSIM, h = Hc)
  to_m <- function(r) if (is_q) -log(1 - r) else r
  ptm <- to_m(fc$rates); colnames(ptm) <- fut
  point <- annuity_e65(make_lookup(win, fut, ptm))
  liab <- numeric(NSIM)
  for (s in 1:NSIM) {
    ms <- to_m(sm$rates[, , s]); colnames(ms) <- fut
    liab[s] <- annuity_e65(make_lookup(win, fut, ms))["annuity"]
  }
  BEL <- mean(liab); SCR <- quantile(liab, 0.995) - BEL
  data.frame(Model = model_name, Window = names(windows)[sapply(windows, identical, win)],
             e65 = round(point["e65"], 2), annuity = round(point["annuity"], 3),
             BEL = round(BEL, 3), SCR = round(SCR, 3),
             SCR_pct = round(100 * SCR / BEL, 1), row.names = NULL)
}

grid <- do.call(rbind, lapply(c("LC", "CBD", "RH"), function(m)
          do.call(rbind, lapply(windows, function(w) run_cell(m, w)))))
print(grid)
write.csv(grid, "figures/capital_grid_R.csv", row.names = FALSE)  # use this to update Table 5.2

grid$Model  <- factor(grid$Model,  levels = c("LC", "CBD", "RH"))
grid$Window <- factor(grid$Window, levels = c("Long", "Recent", "PreSlow"))

# ---- Figure 5.2: best-estimate annuity, grouped bars ----
fig52 <- ggplot(grid, aes(Model, annuity, fill = Window)) +
  geom_col(position = position_dodge(0.72), width = 0.62) +
  coord_cartesian(ylim = c(17, 20)) +
  scale_fill_manual(values = c(NAVY, TEAL, ROSE), name = "window") +
  labs(title = "Best-estimate annuity value by model \u00d7 window",
       x = "Model", y = expression(best*"-"*estimate~annuity~~ddot(a)[65])) +
  theme_diss() + theme(legend.position = "top")
ggsave("figures/fig5_2_annuity_grid.png", fig52, width = 6.4, height = 3.9, dpi = 200)

# ---- Figure 5.3: SCR % of BEL heatmap ----
fig53 <- ggplot(grid, aes(Window, Model, fill = SCR_pct)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.1f", SCR_pct)), colour = "white",
            fontface = "bold", size = 4.6) +
  scale_fill_viridis_c(name = "SCR % of BEL") +
  scale_y_discrete(limits = rev(c("LC", "CBD", "RH"))) +
  labs(title = "Longevity capital: SCR as % of BEL\n(male annuity at 65, by model \u00d7 window)",
       x = "Calibration window", y = "Model") +
  coord_cartesian(expand = FALSE) + theme_diss()
ggsave("figures/fig5_3_capital_heatmap.png", fig53, width = 6.2, height = 3.9, dpi = 200)

cat("\nDone. 7 figures in ./figures/ ; capital grid in figures/capital_grid_R.csv\n")
