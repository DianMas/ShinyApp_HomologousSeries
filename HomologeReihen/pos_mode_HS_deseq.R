install.packages("nontarget")

#MicrOtofQ-Wild potato leaves Exp 9-Stub-Ssto-saca-sbulb_pos.csv
pos_mode_ms_data <- read.delim("pos_mode_msdata.csv.txt", stringsAsFactors = FALSE, header = T, sep = ",")

mz <-pos_mode_ms_data$m.z
rt <- pos_mode_ms_data$RT
e <- pos_mode_ms_data$dummy


count_msdata_posmode <- pos_mode_ms_data[, (-1):(-6)]
row.names(count_msdata_posmode) <- pos_mode_ms_data$Bucket_label
#write.csv(cbind("m/z"=mz, "dummy"=dummy, "RT"=rt),
  #        file="MicrOtofQ-Wild potato leaves Exp 9-Stub-Ssto-saca-sbulb_pos.csv",
   #       row.names=FALSE)
library(readr)
library(nontarget)

# (0.2) list of adducts/isotopes - package enviPat ############
data(adducts);
data(isotopes);

## Use MTBLS1582 peaklist


peaklist <- data.frame(mass=mz, intensity=e, rt=rt)


## according to https://github.com/blosloos/nontarget/issues/6
trace(nontarget::homol.search, function() {}) # put in dummy trace
homol <- homol.search(peaklist,
                      isotopes,	elements = c("C", "H", "O", "N", "P", "S"), use_C = T, # 
                      minmz=12, 	maxmz=200,
                      minrt=-4,  maxrt=-0.5, #0
                      ppm=TRUE,
                      mztol=4,  rttol=0.1, # 0.5
                      minlength=3, #5
                      mzfilter=F,
                      spar=.45, 	R2=.98,
                      plotit=FALSE) # do the thing...
untrace(nontarget::homol.search)
saveRDS(homol, "homol_shiny_app_test.rds")
library(dplyr)
#homol_test <- as.data.frame(homol[["Peaks in homologue series"]][["HS IDs"]] > 0)
#homol.test <- homol[["Peaks in homologue series"]][["HS IDs"]][homol_test,]
#homol_test <- as.data.frame(which(homol[["Peaks in homologue series"]][["HS IDs"]] > 0))
#homol_test <- dplyr::filter(homol[["Peaks in homologue series"]][["peak ID"]], "peak ID" %in% homol_test$`which(homol[["Peaks in homologue series"]][["HS IDs"]] > 0)`)

#homol[["Peaks in homologue series"]][["peak ID"]] <- as.numeric(homol[["Peaks in homologue series"]][["peak ID"]])
library(ggplot2)
# (4.2) Plot results #################################
 plothomol(homol,xlim=FALSE,ylim=FALSE,plotlegend=TRUE, plotdefect = F)

datatest <- as.data.frame(homol[[1]])
###################
binded_homol1_count_msdata <- cbind(homol[[1]], count_msdata_posmode)
write.csv(homol[[3]], file = "homol[[3]]_posmode_3er_possInSource.csv")
#write.csv(homol[[5]], file = "homol[[5]].csv")
write.csv(homol[[1]], file = "homol[[1]]_posmode_3er_-_3001.csv")
write.csv(binded_homol1_count_msdata, file = "homol[[1]]_count_posmode_3er_12-200mzd_-4-05RT.csv")

################################
# Im Folgenden Ansätze für differentiel abundance analysis
################################

homol_metabolites <- read.delim("homol_metabolites.txt", stringsAsFactors = FALSE, header = T, sep = "\t")
meta_msdata <- read.delim("metadata_msdata.txt", stringsAsFactors = F)
count_msdata <- read.delim("count_msdata.txt", stringsAsFactors = F)
#library(mixOmics)
library(DESeq2)
library(edgeR)

#countdata <- count_msdata[,-(1:)]
homol_msdata <- dplyr::filter(count_msdata, m.z %in% homol_metabolites$mz)

#colnames(countdata)
table(colnames(homol_msdata[,6:85])==meta_msdata$sample_Nr)
info_msdata <- homol_msdata[,1:5]
homol_msdata <- homol_msdata[,-(1:5)]
homol_msdata <- data.matrix(homol_msdata, rownames.force = NA)

#rownames(homol_msdata) <- homol_msdata[,3]
#homol_msdata <-  as.matrix(homol_msdata[ , -1])
homol_msdata <- as.matrix(homol_msdata)
#homol_msdata <- homol_msdata + 1
#rownames(meta_msdata) <- meta_msdata[,3]
#head(homol_msdata)
#as.matrix(homol_msdata)

#homol_msdata <- as.numeric(unlist(homol_msdata))        # Apply unlist & as.numeric

is.numeric(homol_msdata)
row.names(homol_msdata) <- info_msdata$Bucket.label
#homol_msdata["" < 0]
dds <- DESeqDataSetFromMatrix(countData = homol_msdata,
                              colData = meta_msdata,
                              design= ~ Species * Treatment)
dds$group <- factor(paste0(dds$Species, dds$Treatment))
design(dds) <- ~ group
dds <- DESeq(dds)
#saveRDS(dds, file ="dds.rds")
resultsNames(dds)
res_homo_SbulboPinf_StuberosumPinf <-  results(dds, contrast = c("group", "StuberosumPinf", "SbulboPinf"), cooksCutoff=T)
res_homo_SbulboPinf_StuberosumPinf <- res_homo_SbulboPinf_StuberosumPinf[order(res_homo_SbulboPinf_StuberosumPinf$padj),]
which_Stolo0_UE0 <- which(res_homo_SbulboPinf_StuberosumPinf$padj < 0.01)# 
which_Stolo0_UE0
#test <- filter(res_homo_SbulboPinf_StuberosumPinf, m.z %in% homol_metabolites$mz)


write.table(res_homo_SbulboPinf_StuberosumPinf[1:144, -(3:4)], file = "res_homo_SbulboPinf_StuberosumPinf.txt", sep = "\t", row.names = T, col.names = T)


#remotes::install_github("rickhelmus/patRoon")
