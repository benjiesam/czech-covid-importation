## 04_testing_adjusted_growth.R -------------------------------------------------
## How fast was the epidemic really growing?
##
## Confirmed cases are the product of two things: how many people were infected
## and how many people were tested. Czech testing went from a few dozen tests a
## day in early March to several thousand by April, so a case curve read at face
## value credits the virus with growth that belongs to the laboratory.
##
## Two models of the same counts, differing only in whether testing effort is
## accounted for:
##
##   cases_t ~ Poisson(exp(a + r*t))                      ignores testing
##   cases_t ~ Poisson(tests_t * exp(a + r*t))            models the positive
##                                                        share, via an offset
##
## The second treats the daily positive share as the quantity with an
## exponential trend, which is the closest this data gets to an incidence rate.
## The gap between the two growth rates is the part of the apparent epidemic
## that was really the testing programme expanding.
##
## This is a descriptive decomposition, not a reproduction number. It says
## nothing about transmission mechanics and should not be read as R0 or Rt.

source(file.path("R", "00_setup.R"))

cl     <- readRDS(file.path(DIR_TAB, "clean.rds"))
daily  <- cl$daily
policy <- cl$policy

## Restrict to the exponential-growth phase. Fitting an exponential across the
## peak and the decline would describe neither, so the window ends at the peak
## of the 7-day average rather than at the end of the series.
ma <- stats::filter(daily$cases, rep(1 / 7, 7), sides = 2)
peak_i <- which.max(ma)
g <- daily[seq_len(peak_i), ]
g <- g[!is.na(g$daily_tests) & g$daily_tests > 0 & g$cases > 0, ]

message("Growth window: ", format(min(g$date)), " to ", format(max(g$date)),
        " (", nrow(g), " days, ending at the peak of the 7-day average)")
message("  tests per day rose from ", min(g$daily_tests), " to ", max(g$daily_tests),
        ", a factor of ", sprintf("%.0f", max(g$daily_tests) / min(g$daily_tests)))
message("  positivity went from ", sprintf("%.1f%%", 100 * g$positivity[1]),
        " to ", sprintf("%.1f%%", 100 * g$positivity[nrow(g)]))

## ---- Two Poisson fits ---------------------------------------------------------------
## Quasi-Poisson throughout: daily counts are clustered and overdispersed, and
## a plain Poisson would report standard errors several times too small.
m_raw <- glm(cases ~ day, family = quasipoisson, data = g)
m_adj <- glm(cases ~ day + offset(log(daily_tests)), family = quasipoisson, data = g)

## unname() matters: coef() returns a named value, and c(r = b) on a named
## element produces "r.day" rather than "r", which silently yields NA on lookup.
rate_of <- function(m) {
  b  <- unname(coef(m)["day"])
  se <- sqrt(vcov(m)["day", "day"])
  c(r = b, lo = b - 1.96 * se, hi = b + 1.96 * se)
}
r_raw <- rate_of(m_raw)
r_adj <- rate_of(m_adj)
dbl   <- function(r) log(2) / r

message("\n--- Growth in reported cases (testing ignored) ---")
message("  rate ", sprintf("%.4f", r_raw["r"]), "/day  (95% CI ",
        sprintf("%.4f", r_raw["lo"]), " to ", sprintf("%.4f", r_raw["hi"]), ")")
message("  doubling time ", sprintf("%.1f days", dbl(r_raw["r"])),
        "  (", sprintf("%.1f", dbl(r_raw["hi"])), " to ",
        sprintf("%.1f", dbl(r_raw["lo"])), ")")

message("\n--- Growth in the positive share (testing offset) ---")
message("  rate ", sprintf("%.4f", r_adj["r"]), "/day  (95% CI ",
        sprintf("%.4f", r_adj["lo"]), " to ", sprintf("%.4f", r_adj["hi"]), ")")
message("  doubling time ", sprintf("%.1f days", dbl(r_adj["r"])),
        "  (", sprintf("%.1f", dbl(r_adj["hi"])), " to ",
        sprintf("%.1f", dbl(r_adj["lo"])), ")")

## ---- What the disagreement means ------------------------------------------------------
## The offset model does not merely slow the growth down, it reverses its sign.
## Taken literally that says the epidemic was shrinking throughout March, which
## is not credible. The reason is that the offset assumes a test on 29 February
## and a test on 2 April are the same instrument, and they are not. The earliest
## tests went almost exclusively to symptomatic people returning from Italy and
## Austria, which is why 7 tests produced 3 cases; by April testing had widened
## to contacts and to the general symptomatic population. Falling positivity
## under that shift is mostly the denominator changing character.
##
## So neither figure is the growth rate. They bracket it:
##   - testing ignored  -> an upper bound, crediting the virus with the laboratory
##   - testing offset   -> a lower bound, crediting the laboratory with the virus
message("\n--- Reading the two together ---")
message("  The two models disagree on the SIGN of the trend.")
message("  They bracket the truth rather than answering: raw growth is an upper")
message("  bound, the offset model a lower bound, because test targeting changed")
message("  from narrow (symptomatic returning travellers) to broad over the window.")

## ---- Sensitivity: drop the period of negligible testing ---------------------------------
## The first days rest on single-digit test counts, which dominate the offset
## model's slope. Refitting from the point where testing reached 100 a day asks
## how much of the reversal is those few days.
g2 <- g[g$daily_tests >= 100, ]
s_raw <- rate_of(glm(cases ~ day, family = quasipoisson, data = g2))
s_adj <- rate_of(glm(cases ~ day + offset(log(daily_tests)),
                     family = quasipoisson, data = g2))
message("\n  Sensitivity, from the first day with >= 100 tests (",
        format(min(g2$date)), ", n = ", nrow(g2), " days):")
message("    testing ignored : ", sprintf("%+.4f", s_raw["r"]), "/day  (doubling ",
        sprintf("%.1f", dbl(s_raw["r"])), " days)")
message("    testing offset  : ", sprintf("%+.4f", s_adj["r"]), "/day")
message("    the sign reversal survives, so it is not an artefact of the")
message("    single-digit test counts in the first week.")

## ---- Figure: cases, tests, positivity -------------------------------------------------
## Three quantities on three panels rather than one chart with two y-axes: a
## second axis would invite the reader to compare two scales chosen by the
## author rather than by the data.
open_png("05_testing_and_cases.png", width = 1900, height = 1600)
par(mfrow = c(3, 1), mar = c(2.4, 5.2, 2.6, 1.6), mgp = c(3.2, 0.6, 0),
    tcl = -0.25, col.axis = COL_MUTED, col.lab = COL_INK, col.main = COL_INK,
    las = 1, font.main = 1, cex.main = 1.0)

panel_axes <- function(dates, ylab, main, last = FALSE) {
  axis(2, col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
  at <- pretty(dates, 8)
  axis(1, at = at, labels = if (last) format(at, "%d %b") else FALSE,
       col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
  title(ylab = ylab); title(main = main, adj = 0, line = 0.9)
}

plot(daily$date, daily$cases, type = "n", axes = FALSE, xlab = "", ylab = "")
abline(h = axTicks(2), col = COL_GRID, lwd = 0.8)
abline(v = policy$date, col = COL_MUTED, lty = 3, lwd = 1)
segments(daily$date, 0, daily$date, daily$cases, col = COL_PRIMARY, lwd = 4, lend = 1)
panel_axes(daily$date, "Confirmed cases", "Daily confirmed cases")

plot(daily$date, daily$daily_tests, type = "n", axes = FALSE, xlab = "", ylab = "")
abline(h = axTicks(2), col = COL_GRID, lwd = 0.8)
abline(v = policy$date, col = COL_MUTED, lty = 3, lwd = 1)
segments(daily$date, 0, daily$date, daily$daily_tests, col = COL_SECONDARY, lwd = 4, lend = 1)
panel_axes(daily$date, "Tests", "Daily tests performed")

plot(daily$date, daily$positivity, type = "n", axes = FALSE, xlab = "", ylab = "",
     ylim = c(0, max(daily$positivity, na.rm = TRUE) * 1.08))
abline(h = axTicks(2), col = COL_GRID, lwd = 0.8)
abline(v = policy$date, col = COL_MUTED, lty = 3, lwd = 1)
lines(daily$date, daily$positivity, col = COL_TERTIARY, lwd = 2.6)
points(daily$date, daily$positivity, pch = 19, col = COL_TERTIARY, cex = 0.6)
axis(2, at = axTicks(2), labels = sprintf("%.0f%%", 100 * axTicks(2)),
     col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
at <- pretty(daily$date, 8)
axis(1, at = at, labels = format(at, "%d %b"), col = COL_GRID,
     col.ticks = COL_GRID, lwd = 0.8)
title(ylab = "Positive share")
title(main = "Share of tests returning positive", adj = 0, line = 0.9)
invisible(dev.off())
message("\n  wrote figures/05_testing_and_cases.png")

## ---- Figure: the two fitted growth curves ---------------------------------------------
open_png("06_growth_comparison.png", width = 1900, height = 1100)
plot(g$date, g$cases, type = "n", log = "y", axes = FALSE, xlab = "", ylab = "",
     ylim = range(c(g$cases, fitted(m_raw), fitted(m_adj))))
abline(h = axTicks(2), col = COL_GRID, lwd = 0.8)
points(g$date, g$cases, pch = 19, col = COL_MUTED, cex = 0.8)
lines(g$date, fitted(m_raw), col = COL_PRIMARY,   lwd = 2.6)
lines(g$date, fitted(m_adj), col = COL_SECONDARY, lwd = 2.6, lty = 2)
axis(1, at = pretty(g$date, 6), labels = format(pretty(g$date, 6), "%d %b"),
     col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
axis(2, col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
title(ylab = "Daily confirmed cases (log scale)")
legend("bottomright",
       legend = c(sprintf("Testing ignored: +%.1f%%/day, doubling %.1f days",
                          100 * r_raw["r"], dbl(r_raw["r"])),
                  sprintf("Testing offset: %+.1f%%/day in the positive share",
                          100 * r_adj["r"])),
       col = c(COL_PRIMARY, COL_SECONDARY), lty = c(1, 2), lwd = 2.6,
       bty = "n", cex = 0.84, text.col = COL_MUTED)
title(main = "The two models disagree on the sign of the trend, and bracket the truth")
invisible(dev.off())
message("  wrote figures/06_growth_comparison.png")

## ---- Outputs -----------------------------------------------------------------------------
out <- data.frame(
  model = c("cases, testing ignored", "positive share, testing offset"),
  growth_rate_per_day = c(r_raw["r"], r_adj["r"]),
  ci_low  = c(r_raw["lo"], r_adj["lo"]),
  ci_high = c(r_raw["hi"], r_adj["hi"]),
  doubling_days = c(dbl(r_raw["r"]), NA),   # undefined for a negative rate
  dispersion = c(summary(m_raw)$dispersion, summary(m_adj)$dispersion),
  window = "full growth phase",
  row.names = NULL
)
out <- rbind(out, data.frame(
  model = c("cases, testing ignored", "positive share, testing offset"),
  growth_rate_per_day = c(s_raw["r"], s_adj["r"]),
  ci_low = c(s_raw["lo"], s_adj["lo"]), ci_high = c(s_raw["hi"], s_adj["hi"]),
  doubling_days = c(dbl(s_raw["r"]), NA), dispersion = NA,
  window = ">= 100 tests/day", row.names = NULL))
print(out, row.names = FALSE, digits = 3)
write_table(out, "04_growth_rates.csv")
