## 03_age_structure.R ----------------------------------------------------------
## If the Czech epidemic was seeded by travel and then spread locally, the two
## groups of cases should not look alike. Austria and Italy account for just
## over half of all imported cases in a season when Czechs travel to both to
## ski, so importation should select working-age adults. Local spread, once it
## reaches households and care settings, should reach everyone.
##
## Three ways of asking the same question, because any one of them can be
## misread: the two distributions side by side, a formal test of whether they
## differ at all, and a logistic model of the probability that a case of a given
## age was travel-attributed, which puts the comparison on a scale that can be
## reported as a number.

source(file.path("R", "00_setup.R"))

cl <- readRDS(file.path(DIR_TAB, "clean.rds"))
ll <- cl$ll

imp <- ll$age[ll$imported == 1]
loc <- ll$age[ll$imported == 0]

message("Imported cases : n = ", length(imp),
        ", median age ", median(imp), ", mean ", sprintf("%.1f", mean(imp)),
        ", IQR ", paste(quantile(imp, c(.25, .75)), collapse = "-"))
message("Not attributed : n = ", length(loc),
        ", median age ", median(loc), ", mean ", sprintf("%.1f", mean(loc)),
        ", IQR ", paste(quantile(loc, c(.25, .75)), collapse = "-"))

## ---- Do the distributions differ? ------------------------------------------------
## The two-sample Kolmogorov-Smirnov statistic is the largest gap between the
## empirical distribution functions, so it detects differences in shape and not
## only in location. Ages are integers and therefore tied, which makes the exact
## p-value unavailable; the asymptotic one is used and the tie warning is
## expected rather than a problem.
ks <- suppressWarnings(ks.test(imp, loc))
wt <- suppressWarnings(wilcox.test(imp, loc, conf.int = TRUE))

message("\nKolmogorov-Smirnov D = ", sprintf("%.3f", ks$statistic),
        ", p ", format.pval(ks$p.value, digits = 3, eps = 1e-16))
message("Hodges-Lehmann shift in median age: ",
        sprintf("%.1f years", wt$estimate),
        " (95% CI ", sprintf("%.1f", wt$conf.int[1]), " to ",
        sprintf("%.1f", wt$conf.int[2]), ")")

## ---- Probability a case was travel-attributed, by age -------------------------------
## Age enters as a natural cubic spline so the shape is not forced to be
## monotone, and calendar day is included because importation declined sharply
## over the window and age composition changed at the same time. Without that
## adjustment the age effect would partly be a time effect in disguise.
ll$day <- as.numeric(ll$date - min(ll$date))

age_knots <- quantile(ll$age, c(0.1, 0.35, 0.65, 0.9))
spl <- function(x) splines::ns(x, knots = age_knots[2:3],
                               Boundary.knots = age_knots[c(1, 4)])

m_crude <- glm(imported ~ spl(age), family = binomial, data = ll)
m_adj   <- glm(imported ~ spl(age) + day, family = binomial, data = ll)
m_null  <- glm(imported ~ day, family = binomial, data = ll)

lrt_age <- anova(m_null, m_adj, test = "LRT")
message("\nAge effect after adjusting for calendar day: deviance drop ",
        sprintf("%.1f", lrt_age$Deviance[2]), " on ", lrt_age$Df[2], " df, p ",
        format.pval(lrt_age$`Pr(>Chi)`[2], digits = 3, eps = 1e-16))

## Predicted probability across the age range, holding the date at the point
## where importation and local transmission were balanced.
cross <- readRDS(file.path(DIR_TAB, "crossover.rds"))
age_grid <- seq(min(ll$age), 95, by = 1)
nd <- data.frame(age = age_grid, day = as.numeric(cross$t_star))
pr <- predict(m_adj, newdata = nd, se.fit = TRUE)
p_hat <- plogis(pr$fit); p_lo <- plogis(pr$fit - 1.96 * pr$se.fit)
p_hi  <- plogis(pr$fit + 1.96 * pr$se.fit)

peak_age <- age_grid[which.max(p_hat)]
message("Probability of travel attribution peaks at age ", peak_age,
        " (", sprintf("%.0f%%", 100 * max(p_hat)), " on the handover date)")
message("  by age 75 it has fallen to ",
        sprintf("%.0f%%", 100 * p_hat[age_grid == 75]))

## ---- Figure: the two age distributions --------------------------------------------
brk <- seq(0, 105, by = 5)
h_imp <- hist(imp, breaks = brk, plot = FALSE)
h_loc <- hist(loc, breaks = brk, plot = FALSE)

open_png("03_age_distributions.png", width = 2000, height = 1150)
par(mar = c(4.4, 4.8, 3.4, 1.6))

plot(NA, xlim = c(0, 100), ylim = c(0, max(h_imp$density, h_loc$density) * 1.12),
     axes = FALSE, xlab = "", ylab = "")
grid_behind()
mid <- h_imp$mids
## Paired bars, offset so neither hides the other, with a gap between them.
rect(mid - 2.3, 0, mid - 0.15, h_imp$density, col = COL_PRIMARY,   border = NA)
rect(mid + 0.15, 0, mid + 2.3, h_loc$density, col = COL_SECONDARY, border = NA)

axis(1, at = seq(0, 100, 10), col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
axis(2, col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
title(xlab = "Age (years)", ylab = "Density within group")
legend("topright", legend = c(sprintf("Travel-attributed (n = %d)", length(imp)),
                              sprintf("Not attributed (n = %d)", length(loc))),
       fill = c(COL_PRIMARY, COL_SECONDARY), border = NA, bty = "n",
       cex = 0.85, text.col = COL_MUTED)
title(main = sprintf(
  "Imported cases cluster in working age; local spread reaches the old (median %d vs %d)",
  median(imp), median(loc)))
invisible(dev.off())
message("\n  wrote figures/03_age_distributions.png")

## ---- Figure: modelled probability by age -------------------------------------------
open_png("04_age_probability.png", width = 1900, height = 1050)
plot(NA, xlim = range(age_grid), ylim = c(0, max(p_hi) * 1.08),
     axes = FALSE, xlab = "", ylab = "")
grid_behind()
polygon(c(age_grid, rev(age_grid)), c(p_lo, rev(p_hi)),
        col = adjustcolor(COL_PRIMARY, alpha.f = 0.16), border = NA)
lines(age_grid, p_hat, col = COL_PRIMARY, lwd = 2.6)
points(peak_age, max(p_hat), pch = 19, col = COL_SECONDARY, cex = 1.2)
text(peak_age, max(p_hat), paste0(" peak at age ", peak_age), adj = c(0, -0.8),
     cex = 0.78, col = COL_SECONDARY)
axis(1, at = seq(0, 95, 10), col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
axis(2, at = axTicks(2), labels = sprintf("%.0f%%", 100 * axTicks(2)),
     col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
title(xlab = "Age (years)", ylab = "P(case attributed to travel)")
title(main = "Modelled on the handover date, adjusted for calendar day")
invisible(dev.off())
message("  wrote figures/04_age_probability.png")

## ---- Outputs --------------------------------------------------------------------------
summ <- data.frame(
  group  = c("travel-attributed", "not attributed"),
  n      = c(length(imp), length(loc)),
  mean   = c(mean(imp), mean(loc)),
  median = c(median(imp), median(loc)),
  q25    = c(quantile(imp, .25), quantile(loc, .25)),
  q75    = c(quantile(imp, .75), quantile(loc, .75)),
  pct_over_65 = c(100 * mean(imp >= 65), 100 * mean(loc >= 65)),
  pct_under_20 = c(100 * mean(imp < 20), 100 * mean(loc < 20)),
  row.names = NULL
)
print(summ, row.names = FALSE, digits = 3)

write_table(summ, "03_age_summary.csv")
write_table(data.frame(age = age_grid, p_imported = p_hat, lo = p_lo, hi = p_hi),
            "03_age_probability_curve.csv")
