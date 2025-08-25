source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

# --------------------------
# User paths (EDIT THESE)
# --------------------------
geno_rds <- "/path/to/geno.rds"
magic_snp <- "/path/to/magic_snp.txt"
pheno_rds <- "path/to/HSC_PCs.rds"

#Read the genotype data from the correct path
geno = readRDS(geno_rds)
myGM = read.table(magic_snp, head = TRUE)

phenotype = readRDS(pheno_rds)
phenotype = as.data.frame(phenotype)
rownames(phenotype) = phenotype[,1]

#Run the GAPIT function, using the GLM model. There are other models and parameters that can be tweaked depending on the phenotype analysis.
myGAPIT_GLM = GAPIT(
    KI = myKI,
    Y = phenotype,
    GD = geno,
    GM = myGM,
    model = c("GLM","MLM","FarmCPU","BLINK"),
    PCA.total = 0,
    file.output = T,
    group.from = 1,
    group.to=616,
    Multiple_analysis=TRUE
  )
