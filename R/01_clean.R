## 01_clean.R ------------------------------------------------------------------
## Assemble three sources into the tables the analysis uses: an individual-level
## line list of confirmed cases, a daily testing series, and region labels.
##
## A note on what `imported_case` means, because everything downstream depends
## on it. The column is 1 for 626 records and empty for the other 6,062, and the
## country-of-exposure field is populated for exactly the same 626. So the flag
## is not missing data in the usual sense: it records whether contact tracing
## attributed the infection to travel. A case with no flag is one that tracing
## did not attribute to travel, which is not the same as a case that was
## certainly acquired locally. Early in an epidemic, with tracing capacity under
## strain, some genuine importations will sit in the unflagged group. That
## biases the imported share downwards, and it biases it downwards most at the
## point where case counts are growing fastest. The estimates in 02 should be
## read as a lower bound on how long importation mattered.

source(file.path("R", "00_setup.R"))

## ---- Line list ----------------------------------------------------------------
ll <- read.csv(file.path(DATA_RAW, "CZ_COVID19_infected_persons.csv"),
               stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

names(ll) <- c("report_date", "age", "gender", "region_code",
               "imported_case", "country_of_exposure")

ll$date     <- as.Date(ll$report_date)
ll$age      <- as.integer(ll$age)
ll$imported <- ifelse(is.na(ll$imported_case) | trimws(ll$imported_case) == "", 0L, 1L)
ll$country_of_exposure[trimws(ll$country_of_exposure) == ""] <- NA_character_
ll$gender   <- factor(trimws(ll$gender), levels = c("Female", "Male"))

stopifnot(!any(is.na(ll$date)), !any(is.na(ll$age)))
stopifnot(all(ll$imported == 1L | is.na(ll$country_of_exposure)))  # flag and country agree

## ---- Region labels --------------------------------------------------------------
reg <- read.csv(file.path(DATA_RAW, "CZ_region_codes.csv"), stringsAsFactors = FALSE,
                encoding = "UTF-8")
names(reg) <- c("region_name", "region_code")
reg$region_name <- trimws(reg$region_name)
## read.csv does not always tag the strings, and an untagged string is rendered
## byte by byte on the graphics device, which is what turns "Plzeň" into "Plze..".
Encoding(reg$region_name) <- "UTF-8"
ll$region <- reg$region_name[match(ll$region_code, reg$region_code)]
stopifnot(!any(is.na(ll$region)))

## ---- Daily counts from the line list ---------------------------------------------
all_days <- seq(min(ll$date), max(ll$date), by = "day")
daily <- data.frame(date = all_days)
daily$cases    <- as.integer(table(factor(ll$date, levels = all_days)))
daily$imported <- as.integer(tapply(ll$imported, factor(ll$date, levels = all_days), sum))
daily$imported[is.na(daily$imported)] <- 0L
daily$local    <- daily$cases - daily$imported
daily$prop_imported <- ifelse(daily$cases > 0, daily$imported / daily$cases, NA_real_)
daily$day <- as.numeric(daily$date - min(daily$date))

## ---- Testing series ---------------------------------------------------------------
## Published separately from the line list and covering a slightly different
## window, so the two are joined rather than assumed to align.
tests <- read.csv(file.path(DATA_RAW, "CZ_COVID19_cumulative_tests_count.csv"),
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
names(tests) <- c("date", "daily_tests", "tests_total")
tests$date <- as.Date(tests$date)

cases_pub <- read.csv(file.path(DATA_RAW, "CZ_COVID19_cumulative_cases_count.csv"),
                      stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
names(cases_pub) <- c("date", "daily_count", "count_total")
cases_pub$date <- as.Date(cases_pub$date)

daily <- merge(daily, tests[, c("date", "daily_tests", "tests_total")],
               by = "date", all.x = TRUE)
daily <- merge(daily, cases_pub[, c("date", "daily_count", "count_total")],
               by = "date", all.x = TRUE)
daily <- daily[order(daily$date), ]
daily$positivity <- daily$cases / daily$daily_tests

## The two sources disagree slightly; say by how much rather than hiding it.
overlap <- !is.na(daily$daily_count)
disc <- sum(daily$cases[overlap]) - sum(daily$daily_count[overlap])

## ---- Known policy dates ------------------------------------------------------------
## Used only as reference lines on figures and never as model inputs.
policy <- data.frame(
  date  = as.Date(c("2020-03-11", "2020-03-12", "2020-03-16")),
  label = c("Schools closed", "State of emergency", "Borders closed, quarantine"),
  stringsAsFactors = FALSE
)

message("Line list     : ", nrow(ll), " cases, ",
        format(min(ll$date)), " to ", format(max(ll$date)))
message("  imported    : ", sum(ll$imported), " (",
        sprintf("%.1f%%", 100 * mean(ll$imported)), ")")
message("  regions     : ", length(unique(ll$region)))
message("  age range   : ", min(ll$age), " to ", max(ll$age),
        " (median ", median(ll$age), ")")
message("Daily series  : ", nrow(daily), " days; testing data on ",
        sum(!is.na(daily$daily_tests)), " of them")
message("  line list vs published daily counts differ by ", disc,
        " cases over the overlap (", sprintf("%.2f%%", 100 * abs(disc) / sum(daily$cases[overlap])),
        "), consistent with later reporting revisions")

top_exp <- sort(table(ll$country_of_exposure), decreasing = TRUE)[1:6]
message("\nTop countries of exposure among imported cases:")
for (i in seq_along(top_exp))
  message("  ", names(top_exp)[i], ": ", top_exp[i],
          sprintf("  (%.1f%% of imported)", 100 * top_exp[i] / sum(ll$imported)))

write_table(daily, "01_daily_series.csv")
saveRDS(list(ll = ll, daily = daily, reg = reg, policy = policy),
        file.path(DIR_TAB, "clean.rds"))
