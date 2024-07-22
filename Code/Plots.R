library(tidyverse)
library(poolr)
library(spectrolab)
library(RColorBrewer)
library(cowplot)
library(gridGraphics)
library(lme4)
library(ape)
library(reshape2)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

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

# Convert a hexadecimal color to RGB, and make it transparent
RGB <- function(hex_color){
  rgb_values <- col2rgb(hex_color)
  rgb_normalized <- as.numeric(rgb_values) / 255
  RGB <- rgb(rgb_normalized[1], rgb_normalized[2], rgb_normalized[3], alpha = 0.5)
  return(RGB)
}

# Indices-GWAS
Ind = readRDS("~/Desktop/GAPIT/Pheno_indices.rds")
Ind_no_taxa <- Ind[,-1]

# Transform the data from wide to long format
Ind_no_taxa_long <- tidyr::gather(Ind_no_taxa, key = "Index", value = "Value")

# Create the boxplot
Fig_1b <- ggplot(Ind_no_taxa_long, aes(x = Index, y = Value)) +
  geom_boxplot() +
  labs(title = "Distribution of Six Indices",
       x = "Index",
       y = "Value")+
  my_theme

cor_matrix_ind <- cor(Ind_no_taxa)
#heatmap(cor_matrix)
# eigenvalues <- eigen(cor_matrix)$values
Meff_Ind <- meff(cor_matrix_ind, method = "nyholt") # 5

#### Multiple Manhattan plots - Indices
phenotype = readRDS("~/Desktop/GAPIT/Pheno_indices.rds")
phenotype = as.data.frame(phenotype)
rownames(phenotype) = phenotype[,1]
myGM = read.table("~/Desktop/Nicotiana/Genomev3.1_MAGIC_QTL/magic_snp.txt", head = TRUE)

setwd("~/Desktop/GAPIT/Indices/")
GMM=GAPIT.Multiple.Manhattan(model_store=c("GLM","MLM","FarmCPU","Blink"),
                             Y.names=c("NDWI","CIre","CCI","ARDSI_Cab","ARDSI_Cw","ARDSI_Cm"),
                             byTraits=TRUE,cutOff=1.83942,plot.type=c("h","s"),GM=myGM)

GMM_CCI=GAPIT.Multiple.Manhattan(model_store=c("GLM","MLM","FarmCPU","Blink"),
                             Y.names="CCI",cutOff=1.83942,plot.type=c("h","s"),GM=myGM)
GAPIT.Circle.Manhattan.Plot(GMM_CCI$multip_mapP,plot.type=c("c","q"),xz=GMM$xz,threshold=1.83942)

## 1e-5
Indices <- read_delim("~/Desktop/GAPIT/Indices_multimodel/filtered_p_e5.txt",skip=1, col_names = FALSE, delim = ",")
Indices = Indices[,c(1:6,8,9)]
colnames(Indices) = c("Wavelength","SNP","Chr","Pos","p_value","MAF","FDR_Adjusted_p_value","effect")

pdf("~/Desktop/P2/Fig/Fig_3b.pdf", width = 7, height = 6)
ggplot(Indices, aes(as.factor(Wavelength), SNP, fill=FDR_Adjusted_p_value)) + 
  geom_tile(width=0.4) +
  scale_fill_gradient(low = "#67000d", high="#fee0d2") +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y = "Chromesome Position",
       x = "Indices",
       fill = "FDR_p-value")
dev.off()

Sig_ind <- filter(Indices, FDR_Adjusted_p_value < 0.05/Meff_Ind) # None

# SW-GWAS
Spec = readRDS("~/Desktop/P2/Pheno_spec_ln.rds")
ml <- read.table("~/Desktop/Nicotiana/Genomev3.1_MAGIC_QTL/cyto.tsv",header = TRUE)

### LMM
Spec_new <- Spec[,-c(3:52)] %>%
  left_join(.,ml,by=c("Taxa" ="Line"))
Spec_new$Line <- as.factor(sub(".*_", "", Spec_new$Taxa))
Spec_new$Cytoplasm = as.factor(Spec_new$Cytoplasm)
Spec_new$Leaf_Num = as.factor(Spec_new$Leaf_Num)

# Initialize a list to store the variance components for each wavelength
variance_components <- data.frame()

# Loop over the wavelengths from the 3rd to the 2103rd column
for (i in 3:2103) {
  # Fit the model for each wavelength
  model <- lmer(Spec_new[,i] ~ (1|Cytoplasm) + (1|Leaf_Num) + (1|Line), data = Spec_new)
  
  # Extract the variance components
  vc <- as.data.frame(VarCorr(model))

  # Calculate the total variance and the percentage
  total_variance <- sum(vc$vcov)
  vc$Percentage <- (vc$vcov / total_variance) * 100
  
  # Store the results in the list
  variance_components <- rbind(variance_components, t(vc$Percentage))
}

colnames(variance_components) <- c("Replicate", "Maternal line", "Leaf_num", "Residual")

wavelengths <- seq(400, 2500, length.out = 2101)
variance_components$Wavelength <- wavelengths

# Reshape the data to a long format
long_variance_components <- melt(variance_components, id.vars = 'Wavelength', variable.name = 'Component', value.name = 'Percentage')

# Create the stacked bar plot
Fig_vc <- ggplot(long_variance_components, aes(x = factor(Wavelength), y = Percentage, fill = Component)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(breaks = seq(500, 2500, by = 500)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("#56B4E9","#D55E00","#009E73","grey")) + 
  labs(x = "Wavelength (nm)", y = "Percentage of Variance (%)", fill = "Component") +
  my_theme

ggsave("~/Desktop/P2/Fig/Fig_vc.pdf",Fig_vc,width = 15, height = 5)

# Remove the "Taxa" column before calculating distances
Spec_no_taxa <- Spec[,-c(1:52)]

CV_Spec = t(as.data.frame(apply(Spec_no_taxa,2,sd) / apply(Spec_no_taxa,2,mean)))

plot_quantile(as_spectra(Spec_no_taxa), col= RGB("#D55E00"), main="All RILs in this study",total_prob = 1, border=FALSE,xlab="Wavelength (nm)", ylab = "Reflectance | CV",
              ylim=c(0,0.6), xaxs = "i")
plot(mean(as_spectra(Spec_no_taxa)), col= RGB("#D55E00"),lty=2,lwd=2,add=TRUE)
plot(as_spectra(CV_Spec), col="#000000",lty=3,lwd=2, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(Spec_no_taxa), col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(1700,0.6,legend=c("Range of reflectance","Mean of reflectance","CV"),
       col=c(NA,"#D55E00","#000000"), bg="transparent", border="transparent",
       lty=c(NA,2,3), lwd = c(NA,2,2),fill=c(RGB("#D55E00"),NA,NA),
       cex=1,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=0.5)
Fig_1a <- recordPlot()

cor_matrix_spec <- cor(Spec_no_taxa)
#heatmap(cor_matrix)
# eigenvalues <- eigen(cor_matrix)$values
Meff_spec <- meff(cor_matrix_spec, method = "nyholt") #1016

Spec_leafnum <- read_delim("~/Desktop/GAPIT/spec_blink_filtered_p_e5.txt",skip=1, col_names = FALSE, delim = ",")
Spec_leafnum = Spec_leafnum[,c(1:6,8,9)]
colnames(Spec_leafnum) = c("Wavelength","SNP","Chr","Pos","p_value","MAF","FDR_Adjusted_p_value","effect")

Spec_leafnum <- read_delim("~/Desktop/GAPIT/SW_MLM/filtered_p_e5.txt",col_names = TRUE, delim = " ")
colnames(Spec_leafnum) = c("Wavelength","SNP","Chr","Pos","p_value","FDR_Adjusted_p_value")


pdf("~/Desktop/P2/Fig/Fig_4b.pdf", width = 14, height = 13)
pdf("~/Desktop/P2/Fig/Sfig_mlm.pdf", width = 10, height = 5)
ggplot((Spec_leafnum), aes(as.numeric(Wavelength), SNP, color=FDR_Adjusted_p_value)) + 
  geom_rect(aes(xmin=400, xmax=700,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=800, xmax=1300,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=1550, xmax=1800,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=2000, xmax=2400,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_vline(xintercept =c(1000,1800),lty=3)+
  geom_point(shape=15,size=5)+
  xlim(400,2500) +
  scale_color_gradient(low = "#67000d", high="#fee0d2") +
  my_theme+
  labs(y = "Chromesome Position",
       x = "Wavelength",
       colour = "FDR_p-value")
dev.off()


Sig_spec_leafnum <- filter(Spec_leafnum,FDR_Adjusted_p_value < 0.05/Meff_spec)

# MLM
spec_sw_mlm <- read.table("~/Desktop/GAPIT/SW_MLM/filtered_p_e5.txt",header = TRUE)
spec_sw_mlm[which.min(spec_sw_mlm$P.value),]

# HSC-PA
results = readRDS("~/Desktop/GAPIT/HSC_new/results.rds")
segments = readRDS("~/Desktop/GAPIT/HSC_new/segments.rds")
segments$wavelength <- as.numeric(segments$wavelength)

# Calculate the range for each segment
segments_ranges <- segments %>%
  arrange(segment, wavelength) %>% # arrange by segment and wavelength
  group_by(segment) %>%
  mutate(diff_wavelength = c(0, diff(wavelength)), # compute the difference between consecutive wavelengths
         group = cumsum(diff_wavelength > 1)) %>% # create a group for each continuous range
  group_by(segment, group) %>%
  summarise(min_wavelength = min(wavelength),
            max_wavelength = max(wavelength)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(Range = ifelse(min_wavelength == max_wavelength, 
                        as.character(min_wavelength), 
                        paste0(min_wavelength, "-", max_wavelength))) %>%
  group_by(segment) %>%
  summarise(Range = paste(unique(Range), collapse=", "))

write.csv(segments_ranges,"~/Desktop/P2/segments_range.csv")

segments_reorganized <- segments %>%
  spread(key = level, value = segment)

# Convert the dataframe to long format
segments_long <- segments_reorganized %>%
  gather(key = "level", value = "segment", -wavelength)

# Convert the level to a factor
segments_long$lwavelength <- as.numeric(segments_long$level)
segments_long$level <- as.factor(segments_long$level)
segments_long$segment <- as.factor(segments_long$segment)

# Create a custom color palette
n_colors <- length(unique(segments_long$segment))
color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(n_colors)
# Create a named vector of colors
named_color_palette <- setNames(color_palette, sort(unique(segments_long$segment)))

# Identify the continuous ranges within each segment
segments_long <- segments_long %>%
  group_by(level, segment) %>%
  mutate(range = cumsum(c(0, diff(as.numeric(wavelength)) > 1)))

# Calculate the middle wavelength for each range
segments_long <- segments_long %>%
  group_by(level, segment, range) %>%
  mutate(middle = round(median(as.numeric(wavelength))))

# Identify segments that span the entire range
segments_spanning <- segments_long %>%
  group_by(level, segment) %>%
  summarize(span = max(as.numeric(wavelength)) - min(as.numeric(wavelength))) %>%
  filter(span == max(span))

# For segments that span the entire range, set the middle wavelength to the median of all wavelengths
segments_long <- segments_long %>%
  left_join(segments_spanning, by = c("level", "segment")) %>%
  mutate(middle = ifelse(!is.na(span), median(as.numeric(wavelength)), middle)) 

# Create the plot
Fig_1c <- ggplot(segments_long, aes(x = as.numeric(wavelength), y = level, fill = segment)) +
  geom_tile(color = NA) +
  #ggrepel::geom_text_repel(data = subset(segments_long, as.numeric(wavelength) == middle), aes(label = segment), size = 5,max.overlaps = Inf) +
  scale_fill_manual(values = named_color_palette,na.value = "transparent") +
  scale_x_continuous(expand = c(0, 0))+
  my_theme +
  labs(title = "Results of HSC-PA",
       x = "Wavelength", y = "Level", fill = "Segment") +
  theme(axis.text.x = element_text(hjust = 1),
        text=element_text(size=20),
        legend.text = element_text(size=15))

# Filter the segments with N_PCs == 1
segments_to_read <- filter(results, N_PCs == 1)$segment

# Plots
final_used = filter(segments,segment %in% segments_to_read)%>%
  mutate(wavelength = as.numeric(wavelength))
example_spectra = apply(Spec_no_taxa,2,mean)
example_spectra_df <- data.frame(wavelength = as.numeric(names(example_spectra)), 
                                 example_spectra = example_spectra)
# Merge the data frames
final_used <- merge(final_used, example_spectra_df, by = "wavelength")

# Create the plot
Fig_1d <- ggplot(final_used, aes(x = wavelength, y = example_spectra, color = as.factor(segment))) +
  geom_rect(aes(xmin=400, xmax=700,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=800, xmax=1300,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=1550, xmax=1800,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=2000, xmax=2400,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_vline(xintercept =c(1000,1800),lty=3)+
  geom_point(size=5) +
  scale_color_manual(values = named_color_palette, name = "Segment") +
  scale_x_continuous(expand = c(0, 0))+
  labs(title = "Segments with 1 retained PC",
       x = "Wavelength (nm)", y = "Reflectance") +
  my_theme

cairo_pdf("~/Desktop/P2/Fig/Fig1.pdf",width = 20,height = 10)
plot_grid(Fig_1a,Fig_1b,Fig_1c,Fig_1d,
          nrow=2,labels = c('(a)','(b)','(c)','(d)'))
dev.off()

## Extract and save the first PC of each phenotypic segment
# Function to read the first column of a file
read_first_column <- function(segment) {
  data <- read.table(paste0("Seg", segment, "_pcs.txt"), header = TRUE)
  return(data[,1])
}

setwd("~/Desktop/GAPIT/HSC_new")
# Read and combine the first column of the files
combined_data <- data.frame(lapply(segments_to_read, read_first_column))
names(combined_data) <- paste0("Seg", segments_to_read)

# Add a new column to combined_data with the taxa information
combined_data$Taxa <- Spec$Taxa

# Move the Taxa column to the first position
combined_data <- combined_data[, c("Taxa", setdiff(colnames(combined_data), "Taxa"))]
saveRDS(combined_data,"~/Desktop/GAPIT/HSC_PCs.rds")

# From here one can also read from the saved RDS file for the HSC-PCs as phenotypes
PCs = readRDS("~/Desktop/GAPIT/HSC_PCs.rds")
PCs_no_taxa <- PCs[,-1]

cor_matrix_PCs <- cor(PCs_no_taxa)
#heatmap(cor_matrix)
# eigenvalues <- eigen(cor_matrix)$values
Meff_PCs <- meff(cor_matrix_PCs, method = "nyholt") #13

Seg <- read_delim("~/Desktop/GAPIT/HSC_final/filtered_p_e5.txt",skip=1, col_names = FALSE, delim = ",")
Seg = Seg[,c(1:6,8,9)]
colnames(Seg) = c("Segments","SNP","Chr","Pos","p_value","MAF","FDR_Adjusted_p_value","effect")

Seg_fdr <- filter(Seg, FDR_Adjusted_p_value < 0.05/Meff_PCs)
Seg_t10 <- arrange(Seg, FDR_Adjusted_p_value) 
Seg_t10 <- Seg_t10[1:10,]

pdf("~/Desktop/P2/Fig/Fig_5b.pdf", width = 10, height = 8)
ggplot(Seg, aes(as.factor(Segments), SNP, fill=FDR_Adjusted_p_value)) + 
  geom_tile(width=0.8) +
  scale_fill_gradient(low = "#67000d", high="#fee0d2") +
  my_theme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y = "Chromesome Position",
       x = "Segments",
       fill = "FDR_p-value")
dev.off()

# Convert segment names in Seg dataframe to numeric format
Seg$Segments <- as.numeric(gsub("Seg", "", Seg$Segments))

# Merge Seg dataframe with final_used dataframe
HSC_fig = merge(final_used, Seg, by.x = "segment", by.y = "Segments")

# Create a column to distinguish between the two datasets
Seg$Type <- "Segment"
HSC_fig$Type <- "Wavelength"

combined_data <- rbind(Seg[,c(2,5,9)], HSC_fig[,c(5,8,12)])

pdf("~/Desktop/P2/Fig/Fig_5c.pdf", width = 16, height = 8)
ggplot(HSC_fig, aes(as.numeric(wavelength), SNP, color=FDR_Adjusted_p_value)) + 
  geom_rect(aes(xmin=400, xmax=700,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=800, xmax=1300,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=1550, xmax=1800,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=2000, xmax=2400,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_vline(xintercept =c(1000,1800),lty=3)+
  geom_point(shape=15,size=5)+
  xlim(400,2500) +
  scale_color_gradient(low = "#67000d", high="#fee0d2") +
  my_theme+
  labs(y = "Chromesome Position",
       x = "Wavelength",
       colour = "FDR_p-value")
dev.off()

test = left_join(final_used,Seg,by = "segment")

HSC_27 = sort(as.numeric(filter(final_used,segment == 27)$wavelength)) # Wavelength 444-503

#### Multiple Manhattan plots - HSC
phenotype = readRDS("~/Desktop/GAPIT/HSC_PCs.rds")
phenotype = as.data.frame(phenotype)
rownames(phenotype) = phenotype[,1]
myGM = read.table("~/Desktop/Nicotiana/Genomev3.1_MAGIC_QTL/magic_snp.txt", head = TRUE)

setwd("~/Desktop/GAPIT/HSC_final/")
GMM=GAPIT.Multiple.Manhattan(model_store=c("GLM","MLM","FarmCPU","Blink"),Y.names="Seg27",cutOff=0.05, GM=myGM)
GAPIT.Circle.Manhattan.Plot(GMM$multip_mapP,plot.type=c("c","q"),xz=GMM$xz,threshold=0.05)

GMM=GAPIT.Multiple.Manhattan(model_store=c("GLM","MLM","FarmCPU","Blink"),
                             Y.names=paste0("Seg",segments_to_read),
                             byTraits=TRUE,cutOff=1.83942,plot.type=c("h","s"),GM=myGM)


setwd("~/Desktop/GAPIT/Results")
# Make plots of sig_p-values
Spec_leafnum <- read_delim("Spec_leafnum/filtered_p_e5.txt",skip=1, col_names = FALSE, delim = " ") %>%
  mutate(Model = "Spec_leafnum")
colnames(Spec_leafnum) = c("Wavelength","SNP","Chr","Pos","p_value","FDR_Adjusted_p_value","Model")

Spec <- read_delim("Spec/filtered_p_e5.txt",skip=1, col_names = FALSE, delim = " ") %>%
  mutate(Model = "Spec")
colnames(Spec) = c("Wavelength","SNP","Chr","Pos","p_value","FDR_Adjusted_p_value","Model")

Norm_leafnum <- read_delim("Normspec_leafnum/filtered_p_e5.txt",skip=1, col_names = FALSE, delim = " ")%>%
  mutate(Model = "Norm_leafnum")
colnames(Norm_leafnum) = c("Wavelength","SNP","Chr","Pos","p_value","FDR_Adjusted_p_value","Model")

Norm <- read_delim("Normspec/filtered_p_e5.txt",skip=1, col_names = FALSE, delim = " ")%>%
  mutate(Model = "Norm")
colnames(Norm) = c("Wavelength","SNP","Chr","Pos","p_value","FDR_Adjusted_p_value","Model")

ggplot((Spec_leafnum), aes(as.numeric(Wavelength), SNP, color=p_value)) + 
  geom_rect(aes(xmin=400, xmax=700,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=800, xmax=1300,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=1550, xmax=1800,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=2000, xmax=2400,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_vline(xintercept =c(1000,1800),lty=3)+
  geom_point(shape=15,size=5)+
  xlim(350,2500) +
  scale_color_gradient(low = "black", high="#F2F2F2") +
  theme_classic()+
  theme(text = element_text(size=23))+
  labs(title = "Spec_leafnum",
       y = "Chromesome Position",
       x = "Wavelength",
       colour = "P-value")

# Select only the desired columns and add the "source" column
All = bind_rows(Spec_leafnum,Spec,Norm_leafnum,Norm) %>%
  filter(FDR_Adjusted_p_value < 0.1) %>%
  group_by(Wavelength, SNP, Chr, Pos) %>%
  summarise(Model = paste(Model, collapse = ", ")) %>%
  ungroup()

SNPs = data.frame(SNP = unique(All$SNP)) %>%
  rowwise()%>%
  mutate(Chromosome = strsplit(SNP, "_")[[1]][1],
         pos = as.numeric(strsplit(SNP, "_")[[1]][2]))

SNPs = data.frame(SNP = unique(Seg_t10$SNP)) %>%
  rowwise()%>%
  mutate(Chromosome = strsplit(SNP, "_")[[1]][1],
         pos = as.numeric(strsplit(SNP, "_")[[1]][2]))

annotation = read.table("~/Desktop/Nicotiana/Genomev3.1_MAGIC_QTL/Niatv3.1.description.csv", header=TRUE, sep=",", quote="\"", row.names = NULL, check.names=FALSE, stringsAsFactors = FALSE)

Anno = SNPs %>%
  left_join(annotation, by = "Chromosome") %>%
  filter(Start>=pos-100000 & End<=pos+100000) %>%
  arrange(SNP)

expression = read.table("~/Desktop/Rishav/tissue_abundance.csv", header=TRUE, sep=",", quote="\"", check.names=FALSE, stringsAsFactors = FALSE)
old = colnames(expression)
old[1] = "GeneID"
colnames(expression) = old

leaves = select(expression, GeneID, leaves_untreated,leaves_normalized,leaves_treated)

Final = left_join(Anno, leaves,by = "GeneID") %>%
  filter(leaves_treated !=0 | leaves_untreated != 0)

# Merge the two data frames based on the SNP column
merged_data <- left_join(Seg_t10, Final, by = "SNP") %>%
  select(SNP, Segments,p_value, FDR_Adjusted_p_value,GeneID,Description,Strand)

# For SNPs without candidate genes, replace NAs with blanks
merged_data[is.na(merged_data)] <- ""

write.csv(merged_data,"~/Desktop/P2/candidategene_top10.csv")

ggplot(All, aes(as.numeric(Wavelength), SNP, color=Model)) + 
  geom_rect(aes(xmin=400, xmax=700,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=800, xmax=1300,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=1550, xmax=1800,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_rect(aes(xmin=2000, xmax=2400,ymin= -Inf, ymax =Inf),fill="grey90", color=NA, alpha=0.05)+
  geom_vline(xintercept =c(1000,1800),lty=3)+
  geom_point(shape=15,size=5)+
  xlim(350,2500) +
  #scale_color_gradient(low = "#F2F2F2", high="black") +
  theme_classic()+
  theme(text = element_text(size=23))+
  labs(y = "Chromesome Position",
       x = "Wavelength",
       colour = "Model")


# PCA
PCs <- read_delim("PCA/filtered_fdr_01.txt",skip=1, col_names = FALSE, delim = " ")
colnames(PCs) = c("PC","SNP","Chr","Pos","p_value","FDR_Adjusted_p_value")

ggplot(PCs, aes(as.factor(PC), SNP, color=p_value)) + 
  geom_point(shape=15,size=5)+
  #xlim(1,18) +
  #scale_x_discrete(breaks = 1:18)+
  scale_color_gradient(low = "#F2F2F2", high="black") +
  theme_classic()+
  theme(text = element_text(size=23))+
  labs(title = "PCA_spec_leafnum",
       y = "Chromesome Position",
       x = "PCs",
       colour = "P-value")

find_ranges <- function(vec) {
  vec <- sort(unique(vec))  # Sort and remove duplicates
  start <- vec[1]
  end <- vec[1]
  ranges <- c()
  
  for (i in 2:length(vec)) {
    if (vec[i] == end + 1) {
      end <- vec[i]
    } else {
      if (start == end) {
        ranges <- c(ranges, as.character(start))
      } else {
        ranges <- c(ranges, paste(start, end, sep = "~"))
      }
      start <- vec[i]
      end <- vec[i]
    }
  }
  
  # Handle the last range
  if (start == end) {
    ranges <- c(ranges, as.character(start))
  } else {
    ranges <- c(ranges, paste(start, end, sep = "~"))
  }
  
  return(paste(ranges, collapse = ","))
}

find_ranges(as.numeric(filter(segments,segment==2)$wavelength))
