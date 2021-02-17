#MicrOtofQ-Wild potato leaves Exp 9-Stub-Ssto-saca-sbulb_neg.csv
neg_mode_ms_data <- read.delim("neg_mode_msdata.csv.txt", stringsAsFactors = FALSE, header = T, sep = ",")

mz <-neg_mode_ms_data$m.z
rt <- neg_mode_ms_data$RT
e <- neg_mode_ms_data$dummy
count_msdata_negmode <- neg_mode_ms_data[, (-1):(-6)] #
row.names(count_msdata_negmode) <- neg_mode_ms_data$Bucket_label

#write.csv(cbind("m/z"=mz, "dummy"=dummy, "RT"=rt),
    #      file="MicrOtofQ-Wild potato leaves Exp 9-Stub-Ssto-saca-sbulb_neg.csv",
     #     row.names=FALSE)
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
                      isotopes,	elements = c("C", "H", "O", "N", "P", "S"), use_C = T, #F?
                      minmz=12, 	maxmz=200, #200?
                      minrt=-0.5,  maxrt=0.5, #0.1 # max=4 ...
                      ppm=TRUE,
                      mztol=4,  rttol=0.1, #0.1 = 6 Sekunden
                      minlength=3,
                      mzfilter=F,
                      spar=.45, 	R2=.98,#98
                      plotit=FALSE) # do the thing...
untrace(nontarget::homol.search)


# (4.2) Plot results #################################
plothomol(homol,xlim=FALSE,ylim=FALSE,plotlegend=TRUE, plotdefect = F)

###################
binded_homol1_count_msdata <- cbind(homol[[1]], count_msdata_negmode)
#write.csv(homol[[3]], file = "homol[[3]]_negmode_3er_-.csv")
#write.csv(homol[[1]], file = "homol[[1]]_negmode_3er_-.csv")
write.csv(binded_homol1_count_msdata, file = "homol[[1]]_count_negmode_4er_12-200mzd_possInSource.csv")

################################
# Im Folgenden Ansätze für differentiel abundance analysis
################################

#homol_metabolites <- read.delim("homol_metabolites.txt", stringsAsFactors = FALSE, header = T, sep = "\t")
#meta_msdata <- read.delim("metadata_msdata.txt", stringsAsFactors = F)
#count_msdata <- read.delim("count_msdata.txt", stringsAsFactors = F)
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

##################
adduct<-adduct.search(peaklist,
                      adducts,
                      rttol=0.1,
                      mztol=5, ppm=TRUE,
                      use_adducts=c("M-H", "M+FA-H", "M+Hac-H"), 
                      ion_mode="negative");
# (2.2) plot results #################################
plotadduct(adduct)


##################
library(tibble)
library(dplyr)
#update.packages("dplyr")
library(igraph)
library(rlang)
library(visNetwork) # for visNetwork() and friends
library(networkD3)  # for saveNetwork()

nl <- as_tibble(peaklist)

colnames(homol[["Peaks in homologue series"]][["to ID"]]) 
el <- as_tibble(homol[["Peaks in homologue series"]]) #%>% 
  
el <-  filter(el, el$`to ID` !="0") #%>% 
el <- select(el, c("peak ID", "to ID", "m/z increment", "RT increment"))

g <- el %>% select (c("peak ID", "to ID")) %>% as.matrix %>% graph_from_edgelist
g
## Write to File
#write_graph(g, "test.gml", "gml")
#getwd()

## Some Plotting
data <- toVisNetworkData(g, idToLabel = T)
head(data)
vn <- visNetwork(nodes = data$nodes, 
                 edges = data$edges) #nodes -> peak IDs
vn
#saveNetwork(vn, "vn.html")
