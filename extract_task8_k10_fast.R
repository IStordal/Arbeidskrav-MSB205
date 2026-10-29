library(tidyverse)
library(lubridate)
library(sf)
library(spdep)
library(spatialreg)
library(lmtest)
library(sandwich)
set.seed(2026)
D <- st_read('lucas_full.gpkg', quiet=TRUE)
if (!inherits(D$sdate_dato,'Date')) D$sdate_dato <- as.Date(D$sdate_dato)
D <- D |> mutate(
 year_month=factor(as.character(year_month),levels=sort(unique(as.character(year_month)))),
 age_sold=year(sdate_dato)-yrbuilt,
 solgt_for_ferdig=age_sold<0,
 inntekt_1000=income_pc/1000,
 stories=fct_collapse(stories,two_plus=c('two','two+half','three')),
 across(c(PRE1960PCT,PTRAF,PNPL,PRMP,PTSDF,PM25),~as.numeric(scale(.x)),.names='{.col}_z'))
h1 <- 'log(price) ~ log(TLA) + log(lotsize) + age_sold + solgt_for_ferdig + rooms + beds + baths + halfbaths + garage + stories + wall + year_month'
h2 <- paste(h1,'+ dist_CBD + inntekt_1000 + andel_svart + andel_asiatisk + andel_andre + andel_blandet + andel_arbeidsledig + mwtt + LESSHSPCT + OVER64PCT + LINGISOPCT')
h3 <- paste(h2,'+ PRE1960PCT_z + PTRAF_z + PNPL_z + PRMP_z + PTSDF_z + PM25_z')
f3 <- as.formula(h3)
B <- st_drop_geometry(D)
lm3 <- lm(f3,data=B)
coords <- st_coordinates(D)
nb <- knn2nb(knearneigh(coords,k=10))
lw <- nb2listw(nb,style='W')
durb <- ~ log(TLA)+log(lotsize)+age_sold+inntekt_1000+andel_svart+PM25_z

cat('OLS HC3\n')
print(coeftest(lm3,vcov.=vcovHC(lm3,type='HC3'))[c('log(TLA)','inntekt_1000'),,drop=FALSE])

cat('\nSLX\n')
slx <- lmSLX(f3,data=B,listw=lw,Durbin=durb)
z <- coef(summary(slx)); print(z[grep('TLA|inntekt_1000',rownames(z)),,drop=FALSE])
rm(slx); gc()

getm <- function(m) {s<-summary(m); if(!is.null(s$Coef)) s$Coef else coef(s)}
cat('\nSEM\n')
sem <- errorsarlm(f3,data=B,listw=lw,method='MC')
z <- getm(sem); print(z[grep('TLA|inntekt_1000',rownames(z)),,drop=FALSE]); cat('lambda=',sem$lambda,'\n')
rm(sem); gc()

cat('\nSDEM\n')
sdem <- errorsarlm(f3,data=B,listw=lw,Durbin=durb,method='MC')
z <- getm(sdem); print(z[grep('TLA|inntekt_1000',rownames(z)),,drop=FALSE]); cat('lambda=',sdem$lambda,'\n')
rm(sdem); gc()

cat('\nSAR\n')
sar <- lagsarlm(f3,data=B,listw=lw,method='MC')
z <- getm(sar); print(z[grep('TLA|inntekt_1000',rownames(z)),,drop=FALSE]); cat('rho=',sar$rho,'\n')

cat('\nIMPACTS\n')
W<-as(lw,'CsparseMatrix'); tr<-trW(W,type='MC'); imp<-impacts(sar,tr=tr,R=200); si<-summary(imp,zstats=TRUE,short=TRUE)
M <- cbind(direct=imp$res$direct,zd=si$zmat[,1],indirect=imp$res$indirect,zi=si$zmat[,2],total=imp$res$total,zt=si$zmat[,3])
print(M[grep('TLA|inntekt_1000',rownames(M)),,drop=FALSE])
