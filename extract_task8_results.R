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
lm3 <- lm(f3, data = lucas_full)

coords <- st_coordinates(lucas_full)
nb3 <- knn2nb(knearneigh(coords, k = 3))
nb10 <- knn2nb(knearneigh(coords, k = 10))
lw3 <- nb2listw(nb3, style = "W")
lw10 <- nb2listw(nb10, style = "W")
boligdata <- st_drop_geometry(lucas_full)
durbin <- ~ log(TLA) + log(lotsize) + age_sold + inntekt_1000 + andel_svart + PM25_z

fit_spatial <- function(lw) {
  list(
    slx = lmSLX(f3, data = boligdata, listw = lw, Durbin = durbin),
    sem = errorsarlm(f3, data = boligdata, listw = lw, method = "MC"),
    sdem = errorsarlm(f3, data = boligdata, listw = lw, Durbin = durbin, method = "MC"),
    sar = lagsarlm(f3, data = boligdata, listw = lw, method = "MC")
  )
}
m3 <- fit_spatial(lw3)
m10 <- fit_spatial(lw10)

fmt <- function(x) sprintf("%.8f", x)

sink("task8_results.txt")
cat("=== CLASSES ===\n")
for (nm in names(m10)) cat(nm, ":", paste(class(m10[[nm]]), collapse=","), "\n")

cat("\n=== OLS HC3 SELECTED ===\n")
ct <- coeftest(lm3, vcov. = vcovHC(lm3, type = "HC3"))
for (term in c("log(TLA)", "inntekt_1000")) {
  cat(term, "coef=", fmt(ct[term,1]), " se=", fmt(ct[term,2]), " t=", fmt(ct[term,3]), " p=", format(ct[term,4], scientific=TRUE), "\n", sep="")
}

get_mat <- function(m) {
  s <- summary(m)
  if (!is.null(s$Coef)) return(s$Coef)
  if (!is.null(s$coefficients)) return(s$coefficients)
  if (!is.null(coef(s))) return(coef(s))
  stop("No coefficient matrix")
}

cat("\n=== K10 COEFFICIENT NAMES MATCHING TLA OR INCOME ===\n")
for (nm in names(m10)) {
  mat <- get_mat(m10[[nm]])
  hits <- grep("TLA|inntekt_1000", rownames(mat), value=TRUE)
  cat("--", nm, "--\n")
  for (term in hits) {
    row <- mat[term,,drop=TRUE]
    cat(term, ":", paste(names(row), fmt(as.numeric(row)), collapse=" | "), "\n")
  }
}

cat("\n=== K3 COEFFICIENT NAMES MATCHING TLA OR INCOME ===\n")
for (nm in names(m3)) {
  mat <- get_mat(m3[[nm]])
  hits <- grep("TLA|inntekt_1000", rownames(mat), value=TRUE)
  cat("--", nm, "--\n")
  for (term in hits) {
    row <- mat[term,,drop=TRUE]
    cat(term, ":", paste(names(row), fmt(as.numeric(row)), collapse=" | "), "\n")
  }
}

cat("\n=== SPATIAL PARAMETERS ===\n")
cat("k3 SEM lambda=", fmt(m3$sem$lambda), " SDEM lambda=", fmt(m3$sdem$lambda), " SAR rho=", fmt(m3$sar$rho), "\n", sep="")
cat("k10 SEM lambda=", fmt(m10$sem$lambda), " SDEM lambda=", fmt(m10$sdem$lambda), " SAR rho=", fmt(m10$sar$rho), "\n", sep="")

impact_extract <- function(m, lw, label) {
  W <- as(lw, "CsparseMatrix")
  tr <- trW(W, type = "MC")
  imp <- impacts(m$sar, tr = tr, R = 200)
  s <- summary(imp, zstats = TRUE, short = TRUE)
  bnames <- attr(imp, "bnames")
  clean <- stringr::str_remove(bnames, " dy/dx$| \\(F\\)$")
  cat("\n=== IMPACTS ", label, " ===\n", sep="")
  for (target in c("log(TLA)", "inntekt_1000")) {
    i <- which(clean == target)
    if (!length(i)) { cat(target, " NOT_FOUND names=", paste(clean, collapse=","), "\n"); next }
    cat(target,
        " direct=", fmt(as.numeric(imp$res$direct[i])),
        " z_direct=", fmt(s$zmat[i,1]),
        " indirect=", fmt(as.numeric(imp$res$indirect[i])),
        " z_indirect=", fmt(s$zmat[i,2]),
        " total=", fmt(as.numeric(imp$res$total[i])),
        " z_total=", fmt(s$zmat[i,3]), "\n", sep="")
  }
}
impact_extract(m3, lw3, "K3")
impact_extract(m10, lw10, "K10")
sink()
