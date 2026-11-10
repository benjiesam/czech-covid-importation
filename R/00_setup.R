## 00_setup.R -----------------------------------------------------------------
## Paths, plotting defaults and small helpers shared by the analysis scripts.
##
## Dependencies are deliberately limited to base R plus MASS, which ships with
## every R installation. The analysis uses no tidyverse and no CRAN downloads,
## so it runs unchanged on a fresh machine or a CI runner.

suppressPackageStartupMessages(library(MASS))

## Czech region names carry diacritics (Ústí nad Labem, Plzeň, Hradec Králové).
## Under a C locale R drops them to dots in plot labels, which is a locale
## problem rather than a font one, so set a UTF-8 locale up front and say so if
## none is available. Containers and CI runners frequently default to C.
if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), fixed = TRUE)) {
  for (loc in c("C.UTF-8", "en_US.UTF-8", "en_GB.UTF-8")) {
    ok <- suppressWarnings(Sys.setlocale("LC_CTYPE", loc))
    if (nzchar(ok)) break
  }
  if (!grepl("UTF-8", Sys.getlocale("LC_CTYPE"), fixed = TRUE))
    warning("No UTF-8 locale available; Czech diacritics will not render in figures.",
            call. = FALSE)
}

## Everything is addressed relative to the project root. Run the scripts from
## there (see run_all.R), not from inside R/.
ROOT <- normalizePath(".", mustWork = TRUE)
if (!dir.exists(file.path(ROOT, "R")))
  stop("Run this from the project root, the directory that contains R/ and data/.")

DATA_RAW <- file.path(ROOT, "data", "raw")
DIR_FIG  <- file.path(ROOT, "outputs", "figures")
DIR_TAB  <- file.path(ROOT, "outputs", "tables")

for (d in c(DIR_FIG, DIR_TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## Colour slots. Three categorical hues, checked for colour-vision deficiency
## separation before use; grey carries everything that is not a data series.
COL_PRIMARY   <- "#2a78d6"   # blue
COL_SECONDARY <- "#eb6834"   # orange
COL_TERTIARY  <- "#1baf7a"   # aqua
COL_INK       <- "#0b0b0b"
COL_MUTED     <- "#52514e"
COL_GRID      <- "#dedcd5"
COL_SURFACE   <- "#fcfcfb"

## A plain, low-chrome plotting frame. Grid first so data draws on top of it.
open_png <- function(file, width = 1800, height = 1200, res = 180) {
  ## type = "cairo" so UTF-8 text is handed to the font stack intact.
  png(file.path(DIR_FIG, file), width = width, height = height, res = res,
      bg = COL_SURFACE, type = "cairo")
  par(mar = c(4.2, 4.6, 3.0, 1.4), mgp = c(2.6, 0.6, 0), tcl = -0.25,
      col.axis = COL_MUTED, col.lab = COL_INK, col.main = COL_INK,
      family = "sans", cex.main = 1.05, font.main = 1, las = 1)
  invisible(NULL)
}

grid_behind <- function(h = TRUE, v = FALSE) {
  if (h) abline(h = axTicks(2), col = COL_GRID, lwd = 0.8)
  if (v) abline(v = axTicks(1), col = COL_GRID, lwd = 0.8)
}

## Light axis furniture: no box, thin ticks, muted labels.
axis_pair <- function(xlab = NULL, ylab = NULL, xat = NULL, xlabels = TRUE) {
  if (is.null(xat)) axis(1, col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
  else axis(1, at = xat, labels = xlabels, col = COL_GRID,
            col.ticks = COL_GRID, lwd = 0.8)
  axis(2, col = COL_GRID, col.ticks = COL_GRID, lwd = 0.8)
  title(xlab = xlab, ylab = ylab)
}

write_table <- function(x, file) {
  path <- file.path(DIR_TAB, file)
  write.csv(x, path, row.names = FALSE)
  message("  wrote ", file, " (", nrow(x), " rows)")
  invisible(path)
}
