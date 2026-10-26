# STA304 EDA. Run Rscript analysis.R from the submission folder.
# Input: unchanged copy of the new 57-record CSV and its codebook.
# This script computes statistics, preparation and graphics in R.
options(stringsAsFactors = FALSE, width = 78)
needed <- c("datawizard", "bayestestR", "car", "jsonlite")
if (!all(vapply(needed, requireNamespace, logical(1), quietly = TRUE))) {
  stop("Install the required packages: ", paste(needed, collapse = ", "))
}
dir.create("tables", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
raw <- read.csv("survey_source.csv", na.strings = "", check.names = FALSE,
  fileEncoding = "UTF-8-BOM", colClasses = c(studentID = "character"))
raw$record <- seq_len(nrow(raw))
raw$bad_id <- is.na(raw$studentID) | !nzchar(trimws(raw$studentID))
total <- with(raw, num_breaks * break_duration_min)
raw$bad_time <- !is.na(total) & total >= raw$study_session_min
raw$bad_week <- !is.na(raw$weekly_study_hours) &
  (raw$weekly_study_hours > 168 | raw$weekly_study_hours < 0)
raw$excluded <- raw$bad_id | raw$bad_time
d <- raw[!raw$excluded, ]
d$weekly_study_hours[d$bad_week] <- NA_real_
# No-break respondents have no defined duration of an actual break.
# Do not infer a duration of zero from the contradictory positive answer.
d$break_duration_min[d$num_breaks == 0] <- NA_real_
split_multi <- function(x) {
  lapply(x, function(z) if (is.na(z)) character() else
    strsplit(z, ",\\s*(?![^()]*\\))", perl = TRUE)[[1]])
}
activities <- split_multi(d$break_activities)
has <- function(label) vapply(activities, function(z) label %in% z, logical(1))
labels <- c(study_session_min = "Study session duration (min)",
  num_breaks = "Breaks per study session",
  break_duration_min = "Duration of each break (min)",
  memory_retention = "Next-day retention rating (1-10)",
  weekly_study_hours = "Weekly study time (hours)",
  prior_grade_pct = "Prior grade: band proxy (%)",
  break_helps_efficiency = "Break benefit agreement (1-7)")
ds <- as.data.frame(datawizard::describe_distribution(d[names(labels)],
  centrality = c("mean", "median")))
write.csv(ds, "tables/easystats_summary.csv", row.names = FALSE)
num <- data.frame(Variable = unname(labels[ds$Variable]), n = ds$n,
  Missing = ds$n_Missing, Mean = ds$Mean, Median = ds$Median,
  SD = ds$SD, IQR = ds$IQR,
  Min = vapply(d[ds$Variable], min, numeric(1), na.rm = TRUE),
  Max = vapply(d[ds$Variable], max, numeric(1), na.rm = TRUE))
freq <- function(x, levels, labels = levels) {
  z <- table(factor(x, levels = levels), useNA = "always")
  data.frame(Category = c(labels, "Missing / unavailable"),
    n = as.integer(z), Percent = 100 * as.integer(z) / length(x))
}
year <- freq(d$year_of_study, 1:5, paste("Recorded code", 1:5))
grade <- freq(d$prior_grade_pct, c(55, 65, 75, 85, 95),
  c("50-59%", "60-69%", "70-79%", "80-89%", "90-100%"))
agree <- freq(d$break_helps_efficiency, 1:7,
  c("Strongly disagree", "Disagree", "Slightly disagree", "Neutral",
    "Slightly agree", "Agree", "Strongly agree"))
memory <- freq(d$memory_retention, 1:10, as.character(1:10))
al <- sort(unique(unlist(activities)))
avalues <- c(Intercept = sprintf("%.5f", coef(ols)[1]),
  Slope = sprintf("%.6f", coef(ols)[2]),
  SlopeP = sprintf("%.4f", coeff$p[2]),
  SlopeT = sprintf("%.3f", coeff$t[2]),
  SlopeLo = sprintf("%.5f", ci[2, 1]),
  SlopeHi = sprintf("%.5f", ci[2, 2]),
  RSquared = sprintf("%.6f", sm$r.squared),
  AdjRSquared = sprintf("%.5f", sm$adj.r.squared),
  ResidualSE = sprintf("%.3f", sm$sigma),
  FStat = sprintf("%.4f", unname(sm$fstatistic[1])),
  TenEffect = sprintf("%.4f", 10 * coef(ols)[2]),
  TenLo = sprintf("%.4f", 10 * ci[2, 1]),
  TenHi = sprintf("%.4f", 10 * ci[2, 2]),
  ShapiroP = sprintf("%.5f", sh_ols$p.value),
  NCVP = sprintf("%.5f", ncv$p),
  MaxCook = sprintf("%.2f", max(cooks)),
  MaxHat = sprintf("%.3f", max(hatvalues(ols))),
  RQOneSlope = sprintf("%.5f", coef(rq1)[2]),
  RQOneP = sprintf("%.4f", rq1_sm$coefficients[2, 4]),
  RQOneLo = sprintf("%.5f", rq1_ci[2, 1]),
  RQOneHi = sprintf("%.5f", rq1_ci[2, 2]),
  AnovaF = sprintf("%.3f", avtab[1, "F value"]),
  AnovaP = sprintf("%.4f", avtab[1, "Pr(>F)"]),
  EtaSquared = sprintf("%.3f", avtab[1, "Sum Sq"] / sum(avtab[, "Sum Sq"])))
names(values) <- sub("\\..*$", "", names(values))
writeLines(paste0("\\newcommand{\\", names(values), "}{", values, "}"),
  "tables/report_values.tex")
navy <- "#23425A"; teal <- "#238A8D"; red <- "#B44B45"
pdf("figures/rq1_results.pdf", width = 7.2, height = 2.5)
par(mar = c(4, 4, 1, 1), cex = .85, las = 1)
set.seed(304)
plot(jitter(d$num_breaks, amount = .07),
  jitter(d$break_helps_efficiency, amount = .05), pch = 19, col = navy,
  xlab = "Breaks per study session", ylab = "Benefit agreement (1-7)")
abline(rq1, col = teal, lwd = 2)
dev.off()
pdf("figures/rq1_diagnostics.pdf", width = 7.2, height = 2.5)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), cex = .8, las = 1)
plot(fitted(rq1), residuals(rq1), pch = 19, col = navy,
  xlab = "Fitted agreement", ylab = "Residual", main = "RQ1: residual spread")
abline(h = 0, lty = 2, col = "grey55")
lines(lowess(fitted(rq1), residuals(rq1)), col = teal, lwd = 2)
qqnorm(residuals(rq1), pch = 19, col = navy, main = "RQ1: normal Q-Q")
qqline(residuals(rq1), col = teal, lwd = 2)
dev.off()
pdf("figures/descriptive_distributions.pdf", width = 7.2, height = 2.65)
par(mfrow = c(1, 2), mar = c(4.1, 3.9, 2.2, .7), cex = .82, las = 1)
hist(bt$break_duration_min, breaks = seq(0, 125, 5), col = "#DCE9EE",
  border = "white", main = "A. Break duration (n = 45)",
  xlab = "Duration of each break (minutes)", ylab = "Respondents")
barplot(as.integer(table(factor(d$memory_retention, levels = 1:10))),
  names.arg = 1:10, col = "#DCE9EE", border = NA,
  main = "B. Retention rating (n = 48)",
  xlab = "Next-day retention rating (1-10)", ylab = "Respondents")
dev.off()
pdf("figures/research_results.pdf", width = 7.2, height = 3.65)
par(mar = c(4.2, 4.2, 2.1, .8), cex = .92, las = 1)
grid <- data.frame(break_duration_min = seq(min(bt$break_duration_min),
  max(bt$break_duration_min), length.out = 250))
prediction <- predict(ols, grid, interval = "confidence")
plot(NA, xlim = c(0, 128), ylim = c(1, 10),
  xlab = "Duration of each break (minutes)",
  ylab = "Next-day retention rating (1-10)",
  main = "Break duration and perceived next-day retention (n = 45)")
polygon(c(grid$break_duration_min, rev(grid$break_duration_min)),
  c(prediction[, "lwr"], rev(prediction[, "upr"])),
  col = adjustcolor(teal, alpha.f = .15), border = NA)
lines(grid$break_duration_min, prediction[, "fit"], col = teal, lwd = 2)
points_n <- aggregate(rep(1, nrow(bt)),
  list(Minutes = bt$break_duration_min, Rating = bt$memory_retention), sum)
points(points_n$Minutes, points_n$Rating, pch = 21,
  bg = adjustcolor(navy, alpha.f = .7), col = navy,
  cex = .85 * sqrt(points_n$x))
im <- which.max(cooks)
points(bt$break_duration_min[im], bt$memory_retention[im],
  pch = 21, col = red, bg = red, cex = 1.1)
text(bt$break_duration_min[im], bt$memory_retention[im],
  labels = "Record 25", pos = 2, col = red, cex = .83)
legend("bottomright", c("1 respondent", "3 respondents", "6 respondents"),
  pch = 21, pt.bg = navy, col = navy, pt.cex = .85 * sqrt(c(1, 3, 6)),
  bty = "n", cex = .78, y.intersp = 1.2)
dev.off()
short <- c("Exercise", "Gaming", "Neither", "Both")
pdf("figures/activity_comparison.pdf", width = 7.2, height = 3.2)
par(mar = c(4.2, 4, 2, .8), cex = .9, las = 1)
boxplot(memory_retention ~ activity_profile, data = ad, names = short,
  col = "#E0ECEF", xlab = "Exercise/gaming profile",
  ylab = "Next-day retention (1-10)", ylim = c(1, 10),
  main = "Retention by exclusive activity profile (n = 46)", outline = FALSE)
set.seed(304)
stripchart(memory_retention ~ activity_profile, data = ad,
  method = "jitter", vertical = TRUE, add = TRUE, pch = 19, col = navy)
dev.off()
pdf("figures/model_diagnostics.pdf", width = 7.2, height = 5.6)
par(mfrow = c(2, 2), mar = c(4, 4, 2.1, 1), cex = .85, las = 1)
plot(bt$break_duration_min, residuals(ols), pch = 19, col = navy,
  xlab = "Break duration (minutes)", ylab = "Residual",
  main = "A. Regression: residual pattern")
abline(h = 0, col = "grey55", lty = 2)
lines(lowess(bt$break_duration_min, residuals(ols)), col = teal, lwd = 2)
qqnorm(residuals(ols), pch = 19, col = navy,
  main = "B. Regression: normal Q-Q", ylab = "Residual quantiles")
qqline(residuals(ols), col = teal, lwd = 2)
boxplot(residuals(av) ~ ad$activity_profile, names = short, col = "#E0ECEF",
  xlab = "Exercise/gaming profile", ylab = "ANOVA residual",
  main = "C. ANOVA: residual spread")
abline(h = 0, col = "grey55", lty = 2)
qqnorm(residuals(av), pch = 19, col = navy,
  main = "D. ANOVA: normal Q-Q", ylab = "Residual quantiles")
qqline(residuals(av), col = teal, lwd = 2)
dev.off()
pdf("figures/data_quality.pdf", width = 7.2, height = 3.35)
par(mfrow = c(1, 2), mar = c(4, 4, 2.2, 1), cex = .84, las = 1)
ok <- !is.na(total) & total > 0
plot(raw$study_session_min[ok], total[ok], log = "xy", pch = 19,
  col = ifelse(raw$bad_time[ok], red, navy),
  xlab = "Session minutes (log scale)",
  ylab = "Total break minutes (log scale)", main = "A. Session-time consistency")
abline(0, 1, col = "grey55", lty = 2)
extreme <- which(raw$bad_time & total > 1000)
text(raw$study_session_min[extreme], total[extreme],
  labels = paste("Record", extreme), pos = 2, cex = .8)
wk <- !is.na(raw$weekly_study_hours) & raw$weekly_study_hours > 0
plot(raw$record[wk], raw$weekly_study_hours[wk], log = "y", pch = 19,
  col = ifelse(raw$bad_week[wk], red, navy), xlab = "Source record number",
  ylab = "Weekly study time (hours; log scale)", main = "B. Weekly-time outliers")
abline(h = 168, lty = 2, col = red)
text(18, 168, "168 hours/week", pos = 3, cex = .8, col = red)
dev.off()
result <- list(method_reference = "proposal_reference.txt",
  raw_n = nrow(raw), screened_n = nrow(d),
  break_takers_n = nrow(breakers), regression_n = nrow(bt), anova_n = nrow(ad),
  bad_id_records = raw$record[raw$bad_id],
  timing_records = raw$record[raw$bad_time],
  weekly_outlier_records = raw$record[raw$bad_week],
  numeric_summary = num, coefficients = coeff, model_fit = fit,
  primary_tests = tests, activity_groups = groups, assumptions = checks,
  rq1 = list(coefficients = rq1_coeff, r_squared = rq1_sm$r.squared,
    residual_df = df.residual(rq1), shapiro_p = rq1_shapiro$p.value,
    variance_score_p = rq1_ncv$p),
  anova = list(df1 = avtab[1, "Df"], df2 = avtab[2, "Df"],
    eta_squared = avtab[1, "Sum Sq"] / sum(avtab[, "Sum Sq"])),
  influence = list(max_cook_record = bt$record[im],
    max_cook = max(cooks), max_leverage = max(hatvalues(ols))))
jsonlite::write_json(result, "results.json", auto_unbox = TRUE,
  pretty = TRUE, digits = 15)
capture.output(sessionInfo(), file = "session_info.txt")
capture.output(citation("datawizard"), file = "tables/datawizard_citation.txt")
