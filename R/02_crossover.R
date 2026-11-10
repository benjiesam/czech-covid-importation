## 02_crossover.R --------------------------------------------------------------
## When did the Czech epidemic stop being driven by travel?
##
## Each day's cases split into those contact tracing attributed to travel and
## those it did not. Modelling that split as binomial in time gives a dated
## answer to a question usually answered by eye:
##
##     logit(P(imported on day t)) = b0 + b1 * t
##
## The handover date is where the fitted probability passes one half, that is
## t* = -b0 / b1. Its uncertainty comes from the delta method and, because
## t* is a ratio of coefficients and ratios behave badly, from a parametric
## bootstrap as a second opinion.
##
## The script then asks whether a single straight line in the log-odds is
## adequate, or whether the decline broke at some date. That is a changepoint
## problem, and because the changepoint is not identified by differentiating the
## likelihood it is found by profiling: fit the segmented model at every
## candidate date and keep the best.

source(file.path("R", "00_setup.R"))

cl     <- readRDS(file.path(DIR_TAB, "clean.rds"))
daily  <- cl$daily
policy <- cl$policy
day0   <- min(daily$date)

## Days with no cases carry no information about the split.
d <- daily[daily$cases > 0, ]
message("Modelling ", nrow(d), " days with at least one case, ",
        format(min(d$date)), " to ", format(max(d$date)), ".")

## ---- Single logistic trend -------------------------------------------------------
fit <- glm(cbind(imported, local) ~ day, family = binomial, data = d)
co  <- coef(fit); V <- vcov(fit)

t_star <- -co[1] / co[2]

## Delta method on t* = -b0/b1.
g      <- c(-1 / co[2], co[1] / co[2]^2)
se_t   <- sqrt(as.numeric(t(g) %*% V %*% g))
ci_delta <- t_star + c(-1.96, 1.96) * se_t

## Parametric bootstrap: resample each day's imported count from the fitted
## binomial, refit, and recompute. This keeps the daily denominators fixed and
## propagates sampling error through the whole pipeline rather than linearising.
set.seed(20200316)
B <- 4000L
## Some replicates produce days of all-zero or all-one counts, where the fitted
## probabilities touch the boundary and glm warns. The fit is still usable for
## the ratio we want, so the warning is caught rather than printed 4,000 times.
boot_t <- replicate(B, {
  p  <- fitted(fit)
  yi <- rbinom(nrow(d), d$cases, p)
  f  <- suppressWarnings(
          try(glm(cbind(yi, d$cases - yi) ~ d$day, family = binomial), silent = TRUE))
  if (inherits(f, "try-error")) return(NA_real_)
  cf <- coef(f)
  if (cf[2] >= 0) return(NA_real_)   # a non-decreasing fit has no crossover
  -cf[1] / cf[2]
})
boot_t <- boot_t[is.finite(boot_t)]
ci_boot <- quantile(boot_t, c(0.025, 0.975))

as_date <- function(x) day0 + x
odds_halflife <- log(2) / -co[2]

message("\n--- Single logistic trend ---")
message("  log-odds slope per day : ", sprintf("%+.4f", co[2]),
        "  (odds of importation halve every ", sprintf("%.1f", odds_halflife), " days)")
message("  crossover (p = 0.5)    : ", format(as_date(t_star)))
message("    delta-method 95% CI  : ", format(as_date(ci_delta[1])), " to ",
        format(as_date(ci_delta[2])))
message("    bootstrap 95% CI     : ", format(as_date(ci_boot[1])), " to ",
        format(as_date(ci_boot[2])),
        "  (", length(boot_t), " usable replicates)")

for (i in seq_len(nrow(policy)))
  message("  ", format(policy$date[i]), "  ", policy$label[i],
          sprintf("   [%+.1f days relative to crossover]",
                  as.numeric(policy$date[i] - as_date(t_star))))

## ---- Overdispersion, and what it does to the interval ------------------------------
## Day-to-day scatter exceeds what binomial sampling allows, which is expected:
## cases arrive in clusters, tracing effort varied by day, and the binomial
## treats every case within a day as independent. A quasi-binomial fit leaves
## the point estimate untouched and scales the standard errors by sqrt(phi),
## which is the honest width for the crossover interval.
fit_q <- glm(cbind(imported, local) ~ day, family = quasibinomial, data = d)
phi   <- summary(fit_q)$dispersion
Vq    <- vcov(fit_q)
se_tq <- sqrt(as.numeric(t(g) %*% Vq %*% g))
ci_quasi <- t_star + c(-1.96, 1.96) * se_tq

message("  residual deviance      : ", sprintf("%.1f", fit$deviance),
        " on ", fit$df.residual, " df (ratio ",
        sprintf("%.2f", fit$deviance / fit$df.residual), ")")
message("  quasi-binomial phi     : ", sprintf("%.2f", phi))
message("    quasi 95% CI         : ", format(as_date(ci_quasi[1])), " to ",
        format(as_date(ci_quasi[2])), "   <- the interval to quote")

## ---- Segmented alternative, by profile likelihood ---------------------------------
## Candidate breaks are restricted so that each segment keeps enough days and
## enough cases to estimate a slope.
cand <- d$day[d$day >= d$day[5] & d$day <= d$day[nrow(d) - 4]]

prof <- do.call(rbind, lapply(cand, function(tau) {
  ## Late in the series the fitted share is within rounding of zero, so glm
  ## reports boundary probabilities. That is the data, not a fault in the fit.
  f <- suppressWarnings(try(glm(cbind(imported, local) ~ day + pmax(day - tau, 0),
                                family = binomial, data = d), silent = TRUE))
  if (inherits(f, "try-error")) return(NULL)
  data.frame(tau = tau, date = as_date(tau), deviance = f$deviance)
}))

best  <- prof[which.min(prof$deviance), ]
fit_s <- suppressWarnings(glm(cbind(imported, local) ~ day + pmax(day - best$tau, 0),
                              family = binomial, data = d))

lrt  <- fit$deviance - fit_s$deviance
## The break date is not estimated under the null, so the usual chi-square on
## 1 df does not apply and this p-value is anti-conservative. Reported as a rough
## guide only; the honest reading is the deviance drop itself.
p_rough <- pchisq(lrt, 1, lower.tail = FALSE)

## Profile-likelihood interval for the break: all candidates within 3.84 deviance
## units of the minimum.
in_ci <- prof$date[prof$deviance <= best$deviance + qchisq(0.95, 1)]

message("\n--- Segmented alternative ---")
message("  best break date        : ", format(best$date))
message("  profile 95% interval   : ", format(min(in_ci)), " to ", format(max(in_ci)))
message("  deviance drop vs single line: ", sprintf("%.1f", lrt),
        "  (rough p ", format.pval(p_rough, digits = 2), ", anti-conservative)")
message("  slope before / after   : ",
        sprintf("%+.4f", coef(fit_s)[2]), " / ",
        sprintf("%+.4f", coef(fit_s)[2] + coef(fit_s)[3]), " per day")

## ---- Figure: the handover ----------------------------------------------------------
grid_days <- seq(min(d$day), max(d$day), by = 0.25)
pr <- predict(fit, newdata = data.frame(day = grid_days), se.fit = TRUE)
mu <- plogis(pr$fit)
lo <- plogis(pr$fit - 1.96 * pr$se.fit)
hi <- plogis(pr$fit + 1.96 * pr$se.fit)

open_png("01_imported_share.png", width = 2000, height = 1280)
par(mar = c(4.4, 4.8, 5.8, 1.6))

plot(d$date, d$prop_imported, type = "n", axes = FALSE, xlab = "", ylab = "",
     ylim = c(0, 1), xlim = range(d$date))
grid_behind()

polygon(c(as_date(grid_days), rev(as_date(grid_days))), c(lo, rev(hi)),
        col = adjustcolor(COL_PRIMARY, alpha.f = 0.16), border = NA)

## Policy dates sit behind the data as reference lines. Two of them are one day
## apart, so the labels alternate in height and fan outwards to avoid colliding.
abline(v = policy$date, col = COL_MUTED, lty = 3, lwd = 1)
lab_y   <- c(1.05, 1.15, 1.05)
lab_adj <- list(c(1, 0), c(0, 0), c(0, 0))
for (i in seq_len(nrow(policy)))
  text(policy$date[i], lab_y[i], paste0(" ", policy$label[i], " "),
       adj = lab_adj[[i]], cex = 0.66, col = COL_MUTED, xpd = NA)

abline(h = 0.5, col = COL_GRID, lwd = 1.2)
## Point area in proportion to the day's case count, so a day resting on three
## cases does not read as strongly as one resting on two hundred.
symbols(d$date, d$prop_imported, circles = sqrt(d$cases), inches = 0.085,
        add = TRUE, bg = adjustcolor(COL_PRIMARY, alpha.f = 0.55), fg = COL_SURFACE)
lines(as_date(grid_days), mu, col = COL_PRIMARY, lwd = 2.4)

segments(as_date(t_star), 0, as_date(t_star), 0.5, col = COL_SECONDARY, lwd = 2)
points(as_date(t_star), 0.5, pch = 19, col = COL_SECONDARY, cex = 1.1)
## The bar shows the quasi-binomial interval, the wider and honest one.
arrows(as_date(ci_quasi[1]), 0.5, as_date(ci_quasi[2]), 0.5, code = 3, angle = 90,
       length = 0.04, col = COL_SECONDARY, lwd = 1.8)
text(as_date(t_star), 0.5, paste0("  handover ", format(as_date(t_star), "%d %b")),
     adj = c(0, 1.9), cex = 0.78, col = COL_SECONDARY)

axis(1, at = pretty(d$date, 8), labels = format(pretty(d$date, 8), "%d %b"),
     col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
axis(2, at = seq(0, 1, 0.25), labels = sprintf("%.0f%%", 100 * seq(0, 1, 0.25)),
     col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
title(ylab = "Share of the day's cases attributed to travel")
title(main = "Local transmission overtook importation before the borders closed",
      line = 4.1)
invisible(dev.off())
message("\n  wrote figures/01_imported_share.png")

## ---- Figure: the profile deviance ----------------------------------------------------
open_png("02_changepoint_profile.png", width = 1800, height = 1000)
plot(prof$date, prof$deviance, type = "n", axes = FALSE, xlab = "", ylab = "")
grid_behind()
abline(h = best$deviance + qchisq(0.95, 1), col = COL_MUTED, lty = 2, lwd = 1.1)
lines(prof$date, prof$deviance, col = COL_PRIMARY, lwd = 2.2)
points(best$date, best$deviance, pch = 19, col = COL_SECONDARY, cex = 1.3)
axis(1, at = pretty(prof$date, 8), labels = format(pretty(prof$date, 8), "%d %b"),
     col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
axis(2, col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
title(ylab = "Residual deviance", xlab = "Candidate break date")
text(best$date, best$deviance, paste0("  ", format(best$date, "%d %b")),
     adj = c(0, 1.6), cex = 0.78, col = COL_SECONDARY)
text(max(prof$date), best$deviance + qchisq(0.95, 1), "95% profile cutoff  ",
     adj = c(1, -0.5), cex = 0.72, col = COL_MUTED)
title(main = "Profiling the break date: every candidate refitted, lowest deviance wins")
invisible(dev.off())
message("  wrote figures/02_changepoint_profile.png")

## ---- Outputs --------------------------------------------------------------------------
res <- data.frame(
  quantity = c("crossover date", "quasi-binomial CI low", "quasi-binomial CI high",
               "bootstrap CI low", "bootstrap CI high",
               "odds halving time (days)", "break date",
               "break profile CI low", "break profile CI high",
               "deviance drop, segmented vs linear"),
  value = c(format(as_date(t_star)), format(as_date(ci_quasi[1])),
            format(as_date(ci_quasi[2])), format(as_date(ci_boot[1])),
            format(as_date(ci_boot[2])), sprintf("%.2f", odds_halflife),
            format(best$date), format(min(in_ci)), format(max(in_ci)),
            sprintf("%.1f", lrt)),
  stringsAsFactors = FALSE
)
write_table(res, "02_crossover_estimates.csv")
write_table(prof, "02_changepoint_profile.csv")
saveRDS(list(fit = fit, fit_s = fit_s, t_star = t_star, day0 = day0,
             ci_boot = ci_boot, best = best), file.path(DIR_TAB, "crossover.rds"))
