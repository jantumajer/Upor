library(dplR)
library(bootRes)

chron_VD <- read.rwl ('f:/upor_n/REC_VD/VD_chronologie_AR.rwl', format=c("auto"))
	upoR_VD <- read.rwl ('f:/upor_n/REC_VD/#VD_Nupor_AR.rwl', format=c("auto")) 
chron_TVA <- read.rwl ('f:/upor_n/REC_VD/TVA_chronologie_AR.rwl', format=c("auto"))
	upoR_TVA <- read.rwl ('f:/upor_n/REC_VD/#TVA_upor_AR.rwl', format=c("auto")) 

scPDSI <- read.table("f:/upor_n/dendroclim_vypocet/spei.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
TEMP <- read.table("f:/upor_n/dendroclim_vypocet/temp.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Q <- read.table("f:/upor_n/dendroclim_vypocet/prutok.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
podzemni <- read.table("f:/upor_n/dendroclim_vypocet/podzemni_voda.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

climate <- list (TEMP, scPDSI); timespan <- c(1933,2012); vnames <- c("TEMP","scPDSI")
# Hydroklimatologie
# climate <- list (scPDSI, Q, podzemni); timespan <- c(1966,2012); vnames <- c("scPDSI","Q","podzemni voda")

############################
# Dendroklimatologie - VD

R2 <- dcc(chron_VD[c(11)], climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
R1 <- dcc(chron_VD[c(12)], climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Z <- dcc(chron_VD[c(13)], climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ref_M <- dcc(chron_VD[c(14)], climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ref_U <- dcc(chron(upoR_VD), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ud <- dcc(chron_VD[c(15)], climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)

# Dendroklimatologie - TVA

R2 <- dcc(chron_TVA[c(11)], climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
R1 <- dcc(chron_TVA[c(12)], climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Z <- dcc(chron_TVA[c(13)], climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ref_M <- dcc(chron_TVA[c(14)], climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ref_U <- dcc(chron(upoR_TVA), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ud <- dcc(chron_TVA[c(15)], climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)

# Dendroklimatologie - VLA

R2 <- dcc(chron(rou2_VLA_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
R1 <- dcc(chron(rou1_VLA_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Z <- dcc(chron(zal_VLA_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ref_M <- dcc(chron(mys_VLA_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ref_U <- dcc(chron(upoR_VLA_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ud <- dcc(chron(upoD_VLA_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)

# Dendroklimatologie - TRW

R2 <- dcc(chron(rou2_TRW_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
R1 <- dcc(chron(rou1_TRW_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Z <- dcc(chron(zal_TRW_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ref_M <- dcc(chron(mys_TRW_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ref_U <- dcc(chron(upoR_TRW_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)
Ud <- dcc(chron(upoD_TRW_serie), climate, method = "corr", start = -5, end = 9, timespan = timespan, vnames = vnames, sb = TRUE, boot = TRUE, ci = 0.05)

############################
# Grafy
dcplot(R1, vertical=FALSE)
dcplot(R2, vertical=FALSE)
dcplot(Z, vertical=FALSE)
dcplot(Ud, vertical=FALSE)
dcplot(Ref_M, vertical=FALSE)
dcplot(Ref_U, vertical=FALSE)


############################
#Ulozeni vysledku

R1$ObsNumber <- 1:34; R2$ObsNumber <- 1:34; Z$ObsNumber <- 1:34; Ud$ObsNumber <- 1:34; Ref_M$ObsNumber <- 1:34; Ref_U$ObsNumber <- 1:34;
names(R1)[c(1,2)] <- c("R1_coef","R1_sig"); names(R2)[c(1,2)] <- c("R2_coef","R2_sig"); names(Z)[c(1,2)] <- c("Z_coef","Z_sig"); names(Ud)[c(1,2)] <- c("Ud_coef","Ud_sig"); names(Ref_M)[c(1,2)] <- c("RefM_coef","RefM_sig"); names(Ref_U)[c(1,2)] <- c("RefU_coef","RefU_sig")
merge_clim <- merge(R1[c(1,2,5)], R2[c(1,2,5)], by="ObsNumber"); merge_clim <- merge(merge_clim, Z[c(1,2,5)], by="ObsNumber"); merge_clim <- merge(merge_clim, Ud[c(1,2,5)], by="ObsNumber"); merge_clim <- merge(merge_clim, Ref_M[c(1,2,5)], by="ObsNumber"); merge_clim <- merge(merge_clim, Ref_U[c(1,2,5)], by="ObsNumber")
write.table(merge_clim, "F:/upor_n/REC_VD/VD_TVA_climate.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")

############################
# Strucchange breakpoints

library(strucchange)

breakpoints <- read.table("F:/upor_n/spline/Data/Analysis/zdroj_pro_strucchange.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)


# For-loop
vstup <- NULL; ss <- NULL; bp.seat <- NULL

for (i in (2:25)){
  vstup <- (c(1933:2012))
  vstup <- cbind(vstup, breakpoints[c(i)])
  colnames(vstup) <- c("d", "y")
  ss <- data.frame(vstup)
  print(names(breakpoints[c(i)])); print(bp.seat <- breakpoints(d~y, data=ss, h = 20))
}

####################################################################################
####################################################################################
####################################################################################
# Jednotlive serie - analyza pointeru

#Pripravne kroky
###############

# Vessel density
serie_VD <- read.rwl ('f:/upor_n/REC_VD/VD_serie_AR.rwl', format=c("auto"))
mys_VD_serie <- serie_VD[c(99:126)];rou1_VD_serie <- serie_VD[c(1:24)];rou2_VD_serie <- serie_VD[c(25:51)];upoD_VD_serie <- serie_VD[c(52:74)];zal_VD_serie <- serie_VD[c(75:98)]

# Total vessel area
serie_TVA <- read.rwl ('f:/upor_n/REC_VD/TVA_serie_AR.rwl', format=c("auto"))
mys_TVA_serie <- serie_TVA[c(99:126)];rou1_TVA_serie <- serie_TVA[c(1:24)];rou2_TVA_serie <- serie_TVA[c(25:51)];upoD_TVA_serie <- serie_TVA[c(52:74)];zal_TVA_serie <- serie_TVA[c(75:98)]

# Tree-ring width
serie_TRW <- read.rwl ('f:/upor_n/TRW/trw_komplet_AR_nonR.rwl', format=c("auto"))
mys_TRW_serie <- serie_TRW[c(76:103)];rou1_TRW_serie <- serie_TRW[c(1:24)];rou2_TRW_serie <- serie_TRW[c(25:51)];upoD_TRW_serie <- serie_TRW[c(104:125)];zal_TRW_serie <- serie_TRW[c(52:75)]
upoR_TRW_serie <- read.rwl('f:/upor_n/TRW/SOUHRN_upoR_AR.rwl', format=c("auto"))

# Vessel lumen area
serie_VLA <- read.rwl ('f:/upor_n/VLA/vla_komplet_detrend.rwl', format=c("auto"))
mys_VLA_serie <- serie_VLA[c(13:40)];rou1_VLA_serie <- serie_VLA[c(41:64)];rou2_VLA_serie <- serie_VLA[c(65:91)];upoD_VLA_serie <- serie_VLA[c(116:137)];upoR_VLA_serie <- serie_VLA[c(1:12)];zal_VLA_serie <- serie_VLA[c(92:115)]
upoR_VLA_serie <- read.rwl ('f:/upor_n/VLA/vla_upoR_detrend.rwl', format=c("auto"))

##########################################

p_mys_VD <- pointer(mys_VD_serie, 100*sens1(serie_VD), 50, 2);p_rou1_VD <- pointer(rou1_VD_serie, 100*sens1(serie_VD), 50, 2);p_rou2_VD <- pointer(rou2_VD_serie, 100*sens1(serie_VD), 50, 2);p_upo_VD <- pointer(upo_VD_serie, 100*sens1(serie_VD), 50, 2);p_zal_VD <- pointer(zal_VD_serie, 100*sens1(serie_VD), 50, 2)
names(p_mys_VD)[c(5)] <- c("Mys_VD"); names(p_rou1_VD)[c(5)] <- c("Rou1_VD"); names(p_rou2_VD)[c(5)] <- c("Rou2_VD"); names(p_upo_VD)[c(5)] <- c("UpoD_VD"); names(p_zal_VD)[c(5)] <- c("zal_VD")
p_merge_VD <- merge(p_mys_VD[c(1,5)], p_rou1_VD[c(1,5)], by="Year"); p_merge_VD <- merge(p_merge_VD, p_rou2_VD[c(1,5)], by="Year"); p_merge_VD <- merge(p_merge_VD, p_upo_VD[c(1,5)], by="Year"); p_merge_VD <- merge(p_merge_VD, p_zal_VD[c(1,5)], by="Year")

p_mys_TVA <- pointer(mys_TVA_serie, 100*sens1(serie_TVA), 50, 2);p_rou1_TVA <- pointer(rou1_TVA_serie, 100*sens1(serie_TVA), 50, 2);p_rou2_TVA <- pointer(rou2_TVA_serie, 100*sens1(serie_TVA), 50, 2);p_upo_TVA <- pointer(upo_TVA_serie, 100*sens1(serie_TVA), 50, 2);p_zal_TVA <- pointer(zal_TVA_serie, 100*sens1(serie_TVA), 50, 2)
names(p_mys_TVA)[c(5)] <- c("Mys_TVA"); names(p_rou1_TVA)[c(5)] <- c("Rou1_TVA"); names(p_rou2_TVA)[c(5)] <- c("Rou2_TVA"); names(p_upo_TVA)[c(5)] <- c("UpoD_TVA"); names(p_zal_TVA)[c(5)] <- c("zal_TVA")
p_merge_TVA <- merge(p_mys_TVA[c(1,5)], p_rou1_TVA[c(1,5)], by="Year"); p_merge_TVA <- merge(p_merge_TVA, p_rou2_TVA[c(1,5)], by="Year"); p_merge_TVA <- merge(p_merge_TVA, p_upo_TVA[c(1,5)], by="Year"); p_merge_TVA <- merge(p_merge_TVA, p_zal_TVA[c(1,5)], by="Year")

p_mys_TRW <- pointer(mys_TRW_serie, 100*sens1(serie_TRW), 50, 2);p_rou1_TRW <- pointer(rou1_TRW_serie, 100*sens1(serie_TRW), 50, 2);p_rou2_TRW <- pointer(rou2_TRW_serie, 100*sens1(serie_TRW), 50, 2);p_upoD_TRW <- pointer(upoD_TRW_serie, 100*sens1(serie_TRW), 50, 2);p_upoR_TRW <- pointer(upoR_TRW_serie, 100*sens1(serie_TRW), 50, 2);p_zal_TRW <- pointer(zal_TRW_serie, 100*sens1(serie_TRW), 50, 2)
names(p_mys_TRW)[c(5)] <- c("Mys_TRW"); names(p_rou1_TRW)[c(5)] <- c("Rou1_TRW"); names(p_rou2_TRW)[c(5)] <- c("Rou2_TRW"); names(p_upoD_TRW)[c(5)] <- c("UpoD_TRW"); names(p_upoR_TRW)[c(5)] <- c("UpoR_TRW"); names(p_zal_TRW)[c(5)] <- c("zal_TRW")
p_merge_TRW <- merge(p_mys_TRW[c(1,5)], p_rou1_TRW[c(1,5)], by="Year"); p_merge_TRW <- merge(p_merge_TRW, p_upoR_TRW[c(1,5)], by="Year"); p_merge_TRW <- merge(p_merge_TRW, p_rou2_TRW[c(1,5)], by="Year"); p_merge_TRW <- merge(p_merge_TRW, p_upoD_TRW[c(1,5)], by="Year"); p_merge_TRW <- merge(p_merge_TRW, p_zal_TRW[c(1,5)], by="Year")

p_mys_VLA <- pointer(mys_VLA_serie, 100*sens1(serie_VLA), 50, 2);p_rou1_VLA <- pointer(rou1_VLA_serie, 100*sens1(serie_VLA), 50, 2);p_rou2_VLA <- pointer(rou2_VLA_serie, 100*sens1(serie_VLA), 50, 2);p_upoD_VLA <- pointer(upoD_VLA_serie, 100*sens1(serie_VLA), 50, 2);p_upoR_VLA <- pointer(upoR_VLA_serie, 100*sens1(serie_VLA), 50, 2);p_zal_VLA <- pointer(zal_VLA_serie, 100*sens1(serie_VLA), 50, 2)
names(p_mys_VLA)[c(5)] <- c("Mys_VLA"); names(p_rou1_VLA)[c(5)] <- c("Rou1_VLA"); names(p_rou2_VLA)[c(5)] <- c("Rou2_VLA"); names(p_upoD_VLA)[c(5)] <- c("UpoD_VLA"); names(p_upoR_VLA)[c(5)] <- c("UpoR_VLA"); names(p_zal_VLA)[c(5)] <- c("zal_VLA")
p_merge_VLA <- merge(p_mys_VLA[c(1,5)], p_rou1_VLA[c(1,5)], by="Year"); p_merge_VLA <- merge(p_merge_VLA, p_upoR_VLA[c(1,5)], by="Year"); p_merge_VLA <- merge(p_merge_VLA, p_rou2_VLA[c(1,5)], by="Year"); p_merge_VLA <- merge(p_merge_VLA, p_upoD_VLA[c(1,5)], by="Year"); p_merge_VLA <- merge(p_merge_VLA, p_zal_VLA[c(1,5)], by="Year")

p_merge_ALL <- merge(p_merge_VD, p_merge_TVA, by="Year"); p_merge_ALL <- merge(p_merge_ALL, p_merge_TRW, by="Year"); p_merge_ALL <- merge(p_merge_ALL, p_merge_VLA, by="Year")
write.table(p_merge_ALL, "F:/upor_n/pointery/p_merge_ALL.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=TRUE, na="NA")

####################
#### Dopocet pointeru pro upoR

p_upoR_TRW <- pointer(upoR_TRW_serie, 100*sens1(serie_TRW), 50, 2)
p_upoR_VLA <- pointer(upoR_VLA_serie, 100*sens1(serie_VLA), 50, 2)
p_upoR_TVA <- pointer(upoR_TVA, 100*sens1(serie_TVA), 50, 2)
p_upoR_VD <- pointer(upoR_VD, 100*sens1(serie_VD), 50, 2)


#######################
# Vypocet popisnych statistik

rwi.stats(rou1_TRW_AR_serie[c(111:190),])
rwi.stats(rou2_TRW_serie[c(111:190),])
rwi.stats(zal_TRW_serie[c(111:190),])
rwi.stats(upoD_TRW_serie[c(111:190),])
rwi.stats(mys_TRW_serie[c(111:190),])
rwi.stats(upoR_TRW_serie[c(31:110),])

rwi.stats(rou1_VLA_serie[c(111:190),])
rwi.stats(rou2_VLA_serie[c(111:190),])
rwi.stats(zal_VLA_serie[c(111:190),])
rwi.stats(upoD_VLA_serie[c(111:190),])
rwi.stats(mys_VLA_serie[c(111:190),])
rwi.stats(upoR_VLA_serie[c(31:110),])

rwi.stats(rou1_VD_serie[c(111:190),])
rwi.stats(rou2_VD_serie[c(111:190),])
rwi.stats(zal_VD_serie[c(111:190),])
rwi.stats(upoD_VD_serie[c(111:190),])
rwi.stats(mys_VD_serie[c(111:190),])
rwi.stats(upoR_VD[c(31:110),])

rwi.stats(rou1_TVA_serie[c(111:190),])
rwi.stats(rou2_TVA_serie[c(111:190),])
rwi.stats(zal_TVA_serie[c(111:190),])
rwi.stats(upoD_TVA_serie[c(111:190),])
rwi.stats(mys_TVA_serie[c(111:190),])
rwi.stats(upoR_TVA[c(31:110),])

