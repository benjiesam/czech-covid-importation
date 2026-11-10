## run_all.R -------------------------------------------------------------------
## Reproduce every figure and table in this repository.
##
## Run from the project root:
##     Rscript run_all.R
##
## Requires base R (>= 4.0) only. No CRAN packages are downloaded.

stopifnot(dir.exists("R"), dir.exists(file.path("data", "raw")))

scripts <- c("R/01_clean.R",
             "R/02_crossover.R",
             "R/03_age_structure.R",
             "R/04_testing_adjusted_growth.R",
             "R/05_regional.R")

t0 <- Sys.time()
for (s in scripts) {
  message("\n", strrep("=", 72), "\n", s, "\n", strrep("=", 72))
  source(s, echo = FALSE)
}

message("\n", strrep("=", 72))
message("Done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
        " seconds.")
message("Figures: ", length(list.files("outputs/figures", "[.]png$")),
        " | Tables: ", length(list.files("outputs/tables", "[.]csv$")))
