library(tidyverse)
library(lubridate)
library(sf)
library(spdep)
library(spatialreg)
library(lmtest)
library(sandwich)

set.seed(2026)
lucas_full <- st_read("lucas_full.gpkg", quiet = TRUE)
if (!inherits(lucas_full$sdate_dato, "Date")) lucas_full$sdate_dato <- as.Date(lucas_full$sdate_dato)

lucas_full <- lucas_full |>
  mutate(
    year_month = factor(as.character(year_month), levels = sort(unique(as.character(year_month)))),
    age_sold = year(sdate_dato) - yrbuilt,
    solgt_for_ferdig = age_sold < 0,
    inntekt_1000 = income_pc / 1000,
    stories = fct_collapse(stories, two_plus = c("two", "two+half", "three")),
    across(c(PRE1960PCT, PTRAF, PNPL, PRMP, PTSDF, PM25),
           ~ as.numeric(scale(.x)), .names = "{.col}_z")
  )

hmod1 <- "log(price) ~ log(TLA) + log(lotsize) + age_sold + solgt_for_ferdig + rooms + beds + baths + halfbaths + garage + stories + wall + year_month"
hmod2 <- paste(hmod1, "+ dist_CBD + inntekt_1000 + andel_svart + andel_asiatisk + andel_andre + andel_blandet + andel_arbeidsledig + mwtt + LESSHSPCT + OVER64PCT + LINGISOPCT")
hmod3 <- paste(hmod2, "+ PRE1960PCT_z + PTRAF_z + PNPL_z + PRMP_z + PTSDF_z + PM25_z")
f3 <- as.formula(hmod3)
boligdata <- st_drop_geometry(lucas_full)
lm3 <- lm(f3, data = boligdata)

coords <- st_coordinates(lucas_full)
nb10 <- knn2nb(knearneigh(coords, k = 10))
lw10 <- nb2listw(nb10, style = "W")
durbin <- ~ log(TLA) + log(lotsize) + age_sold + inntekt_1000 + andel_svart + PM25_z

slx <- lmSLX(f3, data = boligdata, listw = lw10, Durbin = durbin)
sem <- errorsarlm(f3, data = boligdata, listw = lw10, method = "MC")
sdem <- errorsarlm(f3, data = boligdata, listw = lw10, Durbin = durbin, method = "MC")
sar <- lagsarlm(f3, data = boligdata, listw = lw10, method = "MC")
mods <- list(SLX = slx, SEM = sem, SDEM = sdem, SAR = sar)

fmt <- function(x) sprintf("%.8f", as.numeric(x))
get_mat <- function(m) {
  s <- summary(m)
  if (!is.null(s$Coef)) return(s$Coef)
  if (!is.null(s$coefficients)) return(s$coefficients)
  z <- coef(s)
  if (is.matrix(z)) return(z)
  stop("No coefficient matrix")
}

sink("task8_results.txt")
cat("=== OLS HC3 SELECTED ===\n")
ct <- coeftest(lm3, vcov. = vcovHC(lm3, type = "HC3"))
for (term in c("log(TLA)", "inntekt_1000")) {
  cat(term, " coef=", fmt(ct[term,1]), " se=", fmt(ct[term,2]),
      " stat=", fmt(ct[term,3]), " p=", format(ct[term,4], scientific=TRUE), "\n", sep="")
}

cat("\n=== K10 SELECTED COEFFICIENTS ===\n")
for (nm in names(mods)) {
  mat <- get_mat(mods[[nm]])
  hits <- grep("TLA|inntekt_1000", rownames(mat), value = TRUE)
  cat("-- ", nm, " --\n", sep="")
  for (term in hits) {
    row <- mat[term,,drop=TRUE]
    cat(term, ": ", paste(names(row), fmt(row), collapse=" | "), "\n", sep="")
  }
}

cat("\n=== K10 SPATIAL PARAMETERS ===\n")
cat("SEM lambda=", fmt(sem$lambda), "\n", sep="")
cat("SDEM lambda=", fmt(sdem$lambda), "\n", sep="")
cat("SAR rho=", fmt(sar$rho), "\n", sep="")

cat("\n=== K10 SAR IMPACTS ===\n")
W <- as(lw10, "CsparseMatrix")
tr <- trW(W, type = "MC")
imp <- impacts(sar, tr = tr, R = 200)
simp <- summary(imp, zstats = TRUE, short = TRUE)
impact_names <- rownames(as.matrix(imp$res$direct))
if (is.null(impact_names)) impact_names <- names(imp$res$direct)
if (is.null(impact_names)) impact_names <- attr(imp, "bnames")
cat("impact_names=", paste(impact_names, collapse=","), "\n", sep="")
for (target in c("log(TLA)", "inntekt_1000")) {
  i <- match(target, impact_names)
  if (is.na(i)) {
    i <- grep(target, impact_names, fixed=TRUE)[1]
  }
  if (is.na(i) || !length(i)) {
    cat(target, " NOT_FOUND\n", sep="")
  } else {
    cat(target,
        " direct=", fmt(imp$res$direct[i]),
        " z_direct=", fmt(simp$zmat[i,1]),
        " indirect=", fmt(imp$res$indirect[i]),
        " z_indirect=", fmt(simp$zmat[i,2]),
        " total=", fmt(imp$res$total[i]),
        " z_total=", fmt(simp$zmat[i,3]), "\n", sep="")
  }
}
sink()
