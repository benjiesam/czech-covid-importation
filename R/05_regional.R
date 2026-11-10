## 05_regional.R ---------------------------------------------------------------
## Czechia is not one epidemic but fourteen, seeded at different times and from
## different directions. Prague has the international airport; the southern and
## western regions border Austria and Germany and sit on the road home from the
## Alps. If importation mattered unevenly, the regional shares should differ by
## more than chance, and the handover from imported to local should not have
## happened everywhere at once.

source(file.path("R", "00_setup.R"))

cl    <- readRDS(file.path(DIR_TAB, "clean.rds"))
ll    <- cl$ll
cross <- readRDS(file.path(DIR_TAB, "crossover.rds"))
day0  <- cross$day0

## ---- Importation share by region --------------------------------------------------
by_reg <- aggregate(cbind(cases = rep(1, nrow(ll)), imported = ll$imported),
                    by = list(region = ll$region), FUN = sum)
by_reg$prop <- by_reg$imported / by_reg$cases

## Exact (Clopper-Pearson) intervals, because several regions have few cases and
## the normal approximation misbehaves near zero.
ci <- t(mapply(function(x, n) {
  b <- binom.test(x, n)$conf.int
  c(lo = b[1], hi = b[2])
}, by_reg$imported, by_reg$cases))
by_reg <- cbind(by_reg, ci)

## Is the spread more than chance? A chi-square test on the region-by-origin
## table answers directly.
tab <- as.matrix(by_reg[, c("imported", "cases")])
tab <- cbind(imported = tab[, 1], local = tab[, 2] - tab[, 1])
cs  <- chisq.test(tab)

message("Imported share by region, ", nrow(by_reg), " regions")
message("  national   : ", sprintf("%.1f%%", 100 * mean(ll$imported)))
message("  range      : ", sprintf("%.1f%%", 100 * min(by_reg$prop)), " (",
        by_reg$region[which.min(by_reg$prop)], ") to ",
        sprintf("%.1f%%", 100 * max(by_reg$prop)), " (",
        by_reg$region[which.max(by_reg$prop)], ")")
message("  chi-square : ", sprintf("%.1f", cs$statistic), " on ", cs$parameter,
        " df, p ", format.pval(cs$p.value, digits = 3, eps = 1e-16))

## ---- Did the handover happen at different times? --------------------------------------
## One logistic fit per region, keeping only those with enough imported cases to
## support a slope. A region with five imported cases cannot date anything.
ll$day <- as.numeric(ll$date - day0)
min_imported <- 20L

reg_fit <- do.call(rbind, lapply(sort(unique(ll$region)), function(rg) {
  s <- ll[ll$region == rg, ]
  if (sum(s$imported) < min_imported) return(NULL)
  f <- suppressWarnings(glm(imported ~ day, family = binomial, data = s))
  b <- coef(f)
  if (is.na(b[2]) || b[2] >= 0) return(NULL)
  t_star <- -b[1] / b[2]
  g  <- c(-1 / b[2], b[1] / b[2]^2)
  se <- sqrt(as.numeric(t(g) %*% vcov(f) %*% g))
  data.frame(region = rg, n = nrow(s), imported = sum(s$imported),
             crossover = day0 + t_star,
             lo = day0 + t_star - 1.96 * se, hi = day0 + t_star + 1.96 * se,
             stringsAsFactors = FALSE)
}))
reg_fit <- reg_fit[order(reg_fit$crossover), ]

message("\nRegional handover dates (regions with >= ", min_imported,
        " imported cases: ", nrow(reg_fit), " of ", nrow(by_reg), ")")
message("  earliest : ", reg_fit$region[1], ", ",
        format(reg_fit$crossover[1]))
message("  latest   : ", reg_fit$region[nrow(reg_fit)], ", ",
        format(reg_fit$crossover[nrow(reg_fit)]))
message("  spread   : ", sprintf("%.1f days", as.numeric(
          max(reg_fit$crossover) - min(reg_fit$crossover))))
border_date <- cl$policy$date[cl$policy$label == "Borders closed, quarantine"]
n_before <- sum(reg_fit$crossover < border_date)
message("  ", n_before, " of ", nrow(reg_fit), " dated regions crossed over before the ",
        "borders closed on ", format(border_date), ";")
message("  the remaining ", nrow(reg_fit) - n_before,
        " crossed after, so the national date of 9 March conceals real variation.")

## ---- Figure: importation share by region -------------------------------------------------
o  <- order(by_reg$prop)
br <- by_reg[o, ]
m  <- nrow(br)

open_png("07_regional_share.png", width = 1800, height = 1150)
par(mar = c(4.4, 12.5, 3.4, 1.6))
plot(NA, xlim = c(0, max(br$hi) * 1.02), ylim = c(m + 0.5, 0.5),
     axes = FALSE, xlab = "", ylab = "")
abline(v = axTicks(1), col = COL_GRID, lwd = 0.8)
abline(v = mean(ll$imported), col = COL_INK, lwd = 1.3)
segments(br$lo, seq_len(m), br$hi, seq_len(m), col = COL_PRIMARY, lwd = 3, lend = 1)
points(br$prop, seq_len(m), pch = 21, bg = COL_SURFACE, col = COL_PRIMARY,
       cex = 0.95, lwd = 1.8)
axis(1, at = axTicks(1), labels = sprintf("%.0f%%", 100 * axTicks(1)),
     col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
axis(2, at = seq_len(m), labels = paste0(br$region, "  (n=", br$cases, ")"),
     las = 1, cex.axis = 0.7, col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
text(mean(ll$imported), 0.6, " national", adj = c(0, 0), cex = 0.7, col = COL_INK)
title(xlab = "Share of the region's cases attributed to travel")
title(main = sprintf("Importation mattered unevenly across regions (chi-square p %s)",
                     format.pval(cs$p.value, digits = 2, eps = 1e-16)))
invisible(dev.off())
message("\n  wrote figures/07_regional_share.png")

## ---- Figure: regional handover dates ---------------------------------------------------
k <- nrow(reg_fit)
open_png("08_regional_crossover.png", width = 1800, height = 950)
par(mar = c(4.4, 12.5, 3.4, 1.6))
plot(NA, xlim = range(c(reg_fit$lo, reg_fit$hi)), ylim = c(k + 0.5, 0.5),
     axes = FALSE, xlab = "", ylab = "")
abline(v = pretty(c(reg_fit$lo, reg_fit$hi), 6), col = COL_GRID, lwd = 0.8)
abline(v = cl$policy$date, col = COL_MUTED, lty = 3, lwd = 1)
segments(reg_fit$lo, seq_len(k), reg_fit$hi, seq_len(k),
         col = COL_PRIMARY, lwd = 3, lend = 1)
points(reg_fit$crossover, seq_len(k), pch = 21, bg = COL_SURFACE,
       col = COL_PRIMARY, cex = 0.95, lwd = 1.8)
at <- pretty(c(reg_fit$lo, reg_fit$hi), 6)
axis(1, at = at, labels = format(at, "%d %b"), col = COL_GRID,
     col.ticks = COL_GRID, lwd = 0.8)
axis(2, at = seq_len(k), labels = reg_fit$region, las = 1, cex.axis = 0.72,
     col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
text(cl$policy$date[3], 0.55, " borders closed", adj = c(0, 0), cex = 0.68,
     col = COL_MUTED)
title(xlab = "Estimated handover date (imported share passes 50%)")
title(main = sprintf(
  "Handover spans %.0f days: %d of %d dated regions crossed before the borders closed",
  as.numeric(max(reg_fit$crossover) - min(reg_fit$crossover)), n_before, nrow(reg_fit)))
invisible(dev.off())
message("  wrote figures/08_regional_crossover.png")

write_table(by_reg[o, ], "05_regional_share.csv")
write_table(reg_fit, "05_regional_crossover.csv")
