#!/usr/bin/env Rscript

## ================================================================
## Spectra adjustment (cleaned, self-contained)
## - Reads Spec (616 plants), Cyto, AZ (779 plants)
## - Computes indices (in-place, same formulas as your code)
## - Builds factors and UT batch covariate ("ba")
## - For each wavelength, fits: y ~ LN + REP + ba + LN:(REP + ba) + DT:CT + MO
##   (terms with <2 levels are skipped automatically to avoid contrasts() errors)
## - Outputs:
##     pheno_raw.txt, UT_raw.txt (optional)
##     pheno.txt (adjusted spectra + adjusted indices)
##     variance partition figure + summary
## NOTE: spectral column names must be exactly "350","351",...,"2500"
## ================================================================
suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(parallel)
  library(tidyverse)
  library(scales)
  library(ggplot2)
})

# Set global parameters for base R plot and theme for ggplot
par(cex.main = 1.5,  # title size
    font.main = 1,   # title font (2 = bold)
    cex.lab = 1.4,   # label size
    font.lab = 1,    # label font (2 = bold)
    cex.axis = 1.4,  # axis text size
    font.axis = 1,
    mar = c(4,4,4,3) + 0.1,
    bg = "transparent")   

my_theme <- theme_classic()+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        text=element_text(size=17),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill="transparent"))

# --------------------------
# User paths (EDIT THESE)
# --------------------------
spec_rds <- "/path/to/Pheno_leafspec.rds"
cyto_tsv <- "/path/to/cyto.tsv"
az_rds   <- "/path/to/R_AZ.rds"
out_dir  <- "path/to/your/output/directory"
fig_dir <- "paht/to/your/figures/directory"

## 1) Read phenotype raw data:
Spec <- readRDS(spec_rds);dim(Spec) #spectral data of 616 plants
Cyto <- fread(cyto_tsv, data.table=FALSE);dim(Cyto) #maternal cytoplasm groups (650 plants)
AZ <- readRDS(az_rds);dim(AZ) #spectral data of 779 plants (including above ones)
AZ$tempTaxa <- gsub("M","L",AZ$Genotype_ID) #replacing M by L in names of individuals
AZ$Taxa <- gsub("-","_",AZ$tempTaxa) #replacing hyphen by underscore
split_Line <- strsplit(AZ$Genotype_ID,"-") #converting names of individuals into relicate and genotype
AZ$genotype_group <- sapply(split_Line,function(x) x[1])
AZ$genotype_number <- sapply(split_Line,function(x) x[2])
az <- subset(AZ,select=c(Plant_ID,Leaf_Num,Batch,Day,Ctime,Taxa,genotype_group,genotype_number))
dim(az) #Plant_ID,Leaf_Num,Batch,Day,Ctime,Taxa,genotype_group,genotype_number for 779 plants
az <- subset(az,genotype_group=="M1"|genotype_group=="M2")
dim(az) # 628 plants
nlevels(factor(az$Taxa)) #620. there are duplicates which are not in Spec, but which ones?
tapply(az$Day,factor(az$Taxa),length) #L1_111,L1_290,L1_76,L1_78,L2_264,L2_182,L2_194,L2_99 have two rows
#remove second occurrence of the above:
az <- filter(az %>% distinct(Taxa,.keep_all = TRUE))
dim(az)
UTonly <- subset(AZ,genotype_group=="UT") #spectra of phytometer plants (Utah wildtype genotype)
dim(UTonly) #these are the 99 UT plants
#merging the files with different information about same plants:
pheno_temp <- merge(Cyto,Spec,by.x="Line",by.y="Taxa")
dim(pheno_temp) #616 plants. There are four plants without genotypes
pheno <- merge(pheno_temp,az,by.x="Line",by.y="Taxa")
dim(pheno) #file with complete information about phenotypes of 616 plants
pheno[c(1:3,614:616),c(1:3,2150:2160)]
nlevels(factor(pheno$Line)) #just checking this is unique

## 2) Calculate indices:
#pheno:
nm865 <- pheno[,865-349+2];nm1614 <- pheno[,1614-349+2]
pheno$NDWI <- (nm865-nm1614)/(nm865+nm1614) 
nm783 <- pheno[,783-349+2];nm704 <- pheno[,704-349+2]
pheno$CIre <- nm783/nm704-1 
nm560 <- pheno[,560-349+2];nm664 <- pheno[,664-349+2]
pheno$CCI <- (nm560-nm664)/(nm560+nm664)
nm750 <- pheno[,750-349+2];nm730 <- pheno[,730-349+2]
nm770 <- pheno[,770-349+2];nm720 <- pheno[,720-349+2]
pheno$ARDSI_Cab <- (nm750-nm730)/(nm770+nm720)
nm1360 <- pheno[,1360-349+2];nm1080 <- pheno[,1080-349+2]
nm1560 <- pheno[,1560-349+2];nm1240 <- pheno[,1240-349+2]
pheno$ARDSI_Cw <- (nm1360-nm1080)/(nm1560+nm1240)
nm2200 <- pheno[,2200-349+2];nm1640 <- pheno[,1640-349+2]
nm2240 <- pheno[,2240-349+2];nm1720 <- pheno[,1720-349+2]
pheno$ARDSI_Cm <- (nm2200-nm1640)/(nm2240+nm1720)
dim(pheno)
pheno <- pheno %>% relocate(NDWI,CIre,,CCI,ARDSI_Cab,ARDSI_Cw,ARDSI_Cm
                            ,.before=Plant_ID)
#UTonly:
nm865 <- UTonly[,865-349];nm1614 <- UTonly[,1614-349]
UTonly$NDWI <- (nm865-nm1614)/(nm865+nm1614) 
nm783 <- UTonly[,783-349];nm704 <- UTonly[,704-349]
UTonly$CIre <- nm783/nm704-1 
nm560 <- UTonly[,560-349];nm664 <- UTonly[,664-349]
UTonly$CCI <- (nm560-nm664)/(nm560+nm664)
nm750 <- UTonly[,750-349];nm730 <- UTonly[,730-349]
nm770 <- UTonly[,770-349];nm720 <- UTonly[,720-349]
UTonly$ARDSI_Cab <- (nm750-nm730)/(nm770+nm720)
nm1360 <- UTonly[,1360-349];nm1080 <- UTonly[,1080-349]
nm1560 <- UTonly[,1560-349];nm1240 <- UTonly[,1240-349]
UTonly$ARDSI_Cw <- (nm1360-nm1080)/(nm1560+nm1240)
nm2200 <- UTonly[,2200-349];nm1640 <- UTonly[,1640-349]
nm2240 <- UTonly[,2240-349];nm1720 <- UTonly[,1720-349]
UTonly$ARDSI_Cm <- (nm2200-nm1640)/(nm2240+nm1720)
dim(UTonly)
UTonly <- UTonly %>% relocate(NDWI,CIre,,CCI,ARDSI_Cab,ARDSI_Cw,ARDSI_Cm
                              ,.before=Genotype_ID)

## 3) Write phenotype raw data # Optional
write.table(pheno,file = file.path(out_dir, "pheno_raw.txt"),quote=F,row.names=F,sep="\t")
write.table(UTonly,file = file.path(out_dir, "UT_raw.txt"),quote=F,row.names=F,sep="\t")

## Make factors and add to dataset pheno (and UTonly):
pheno$BA <- factor(pheno$Batch);nlevels(pheno$BA)
pheno$DT <- factor(pheno$Day);levels(pheno$DT)
pheno$CT <- factor(pheno$Ctime);levels(pheno$CT)
split_Line <- strsplit(pheno$Line,"_")
repl <- sapply(split_Line,function(x) x[1])
plant <- sapply(split_Line,function(x) x[2])
pheno$REP <- factor(repl);nlevels(pheno$REP)
pheno$PL <- factor(plant);nlevels(pheno$PL)
pheno$MO <- factor(pheno$Cytoplasm);nlevels(pheno$MO)
pheno$LN <- factor(pheno$Leaf_Num);nlevels(pheno$LN)
dim(pheno)
UTonly$BA <- factor(UTonly$Batch);nlevels(UTonly$BA) #first batch no UT

## 4) Prepare adjusted phenotypic data for downstream analysis:
## adjusting reflectance values for non-genetic influences using for loop
## (note that this is only done for raw spectral data and indices,
## segment PCs should be calculated from the corrected reflectance data):
## Make factors and add to dataset pheno (and UTonly):
dim(pheno)
names(pheno[2153:2173]) #again checking variables in final columns (last two not needed)
UTonly$BA <- factor(UTonly$Batch);nlevels(UTonly$BA) #first batch no UT
pheno[1,c(1:4,2153:2173)] #column 3 is 350 nm
UTonly[1,1:2] #column 1 is 350 nm
#dimensions of results:
names(pheno[,c(1,3,2153:2159)])
n_wvl <- 2157;n_pl <- length(pheno$PL)
newy <- matrix(data=0,nrow=n_pl,ncol=n_wvl)

terms_vec <- c("LN","REP","ba","LN:REP","LN:ba","DT:CT","MO","Residuals")
ss_list <- vector("list", length = 2159 - 3 + 1)
wl_vec  <- numeric(length(ss_list))
#start of loop:
for(i in 3:2159) {#i <- 3-350+700 to run 700 nm as example
  pheno$y <- pheno[,i]
  #characterize batches by mean values for UT at i-3+350 nm:
  BAmeans <- tapply(UTonly[,i-2],UTonly$BA,mean,na.rm=T)
  UT_BA <- BAmeans[1:34]
  batch1 <- mean(UT_BA)
  utba <- c(batch1,UT_BA) #missing value for first batch replaced by overall mean of UT
  pheno$ba <- utba[pheno$BA] #assigns UT value of each batch
  #extract residuals after fitting environmental terms:
  lmy2 <- lm(terms(pheno$y~LN #influence measurement technique
                   +REP+ba+LN:(REP+ba) #influence of space and m.t.
                   +DT:CT #influence time
                   +MO #Maternal cytoplasm
                   ,keep.order=T)
             ,data=pheno)
  #anova(lmy2)
  new <- residuals(lmy2)+mean(pheno$y)
  newy[,i-2] <- new
  
  a <- anova(lmy2)                              # base R ANOVA (type I SS)
  ss <- setNames(rep(NA_real_, length(terms_vec)), terms_vec)
  
  rn <- rownames(a)
  if ("Residuals" %in% rn) ss["Residuals"] <- a["Residuals","Sum Sq"]
  
  # Fill available terms; note LN:(REP+ba) expands to LN:REP and LN:ba in 'a'
  present <- intersect(terms_vec, rn)
  ss[present] <- a[present, "Sum Sq"]
  
  tot <- sum(ss, na.rm = TRUE)
  props <- ss / tot
  
  ss_list[[i - 2]] <- props
  wl_vec[i - 2] <- i - 3 + 350
  
}

newy[1:10,c(2151:2157)] #environmentally corrected reflectances and indices for each wavelength
## Write corrected phenotypic data for 616 individuals (note genotype 245 is missing):
newy2 <- data.frame(pheno$Line,pheno$REP,pheno$PL,newy)
dim(newy2)
names(newy2)[names(newy2)=="X2152"] <- "NDWI"
names(newy2)[names(newy2)=="X2153"] <- "CIre"
names(newy2)[names(newy2)=="X2154"] <- "CCI"
names(newy2)[names(newy2)=="X2155"] <- "ARDSI_Cab"
names(newy2)[names(newy2)=="X2156"] <- "ARDSI_Cw"
names(newy2)[names(newy2)=="X2157"] <- "ARDSI_Cm"
colnames(newy2)[1:2154] = c("Taxa","REP","PL",seq(350,2500))
write.table(newy2,file=file.path(out_dir,"pheno.txt"),quote=F,row.names=F,sep="\t") #corrected spectra
## Generate adjusted spec file and adjusted indices file
write_rds(newy2[c(1,4:2154)],file = file.path(out_dir,"Pheno_corr.rds"))
write_rds(newy2[c(1,2155:2160)],file = file.path(out_dir,"Indices_corr.rds"))

## 5) Variance partition analyses and plots
var_df <- do.call(rbind, lapply(seq_along(ss_list), function(k) {
  data.frame(wavelength = wl_vec[k],
             term = names(ss_list[[k]]),
             prop = as.numeric(ss_list[[k]]),
             row.names = NULL, check.names = FALSE)
}))
var_df$term <- factor(var_df$term, levels = terms_vec)

## ==== define groups ====
var_df$group <- NA_character_
var_df$group[var_df$term == "LN"] <- "Technique"
var_df$group[var_df$term %in% c("REP", "ba", "LN:REP", "LN:ba")] <- "Space + tech"
var_df$group[var_df$term == "DT:CT"] <- "Time"
var_df$group[var_df$term == "MO"] <- "Maternal cytoplasm"
var_df$group[var_df$term == "Residuals"] <- "Residuals"

# Check no NAs left
stopifnot(!any(is.na(var_df$group)))

## ==== collapse by group per wavelength ====
group_df <- aggregate(prop ~ wavelength + group, data = var_df, sum, na.rm = TRUE)
group_df$group <- factor(group_df$group,
                         levels = c("Space + tech",
                                    "Technique",
                                    "Time",
                                    "Maternal cytoplasm",
                                    "Residuals"))

## ==== plot: 100% stacked area by group ====
Fig_variance_part <- ggplot(group_df, aes(x = wavelength, y = prop, fill = group)) +
  geom_area() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0,0)) +
  scale_x_continuous(breaks = seq(500, 2500, by = 500),expand=c(0,0)) +
  #scale_x_discrete(breaks = seq(500, 2500, by = 500),expand = c(0,0)) +
  scale_fill_manual(
    values = c(
      "Residuals" = "grey70",
      "Space + tech" = "#009E73",
      "Technique" = "#CC79A7",
      "Time" = "#56B4E9",
      "Maternal cytoplasm" = "#D55E00"
    )
  ) +
  labs(x = "Wavelength (nm)",
       y = "Variance explained (% of total)",
       fill = "Effect group") +
  my_theme

ggsave(file.path(fig_dir,"Variance_part.pdf"),Fig_variance_part,width = 15, height = 5)

## ==== optional: mean contribution per group ====
mean_group <- aggregate(prop ~ group, data = group_df, mean, na.rm = TRUE)
ggplot(mean_group, aes(x = group, y = prop, fill = group)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Mean variance explained",
       title = "Average variance contribution per group") +
  theme_bw() +
  theme(legend.position = "none")

