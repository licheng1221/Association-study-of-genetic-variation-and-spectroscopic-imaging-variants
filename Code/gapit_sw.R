source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#Read the genotype data from the correct path
geno = readRDS("../geno.RDS")
myGM = read.table("../magic_snp.txt", head = TRUE)
myKI = read.table("../magic_kinship.txt", head = TRUE)
cv = read.table("../covariate_leafnum.txt", header=TRUE,sep=",", stringsAsFactors=FALSE)
rownames(cv) = cv[,1]

phenotype = readRDS("../HSC_PCs.rds")
phenotype = as.data.frame(phenotype)
rownames(phenotype) = phenotype[,1]

#Run the GAPIT function, using the GLM model. There are other models and parameters that can be tweaked depending on the phenotype analysis.
for (i in 1:(ncol(phenotype)-1)){
  myGAPIT_GLM = GAPIT(
    KI = myKI,
    Y = phenotype,
    GD = geno,
    GM = myGM,
    CV = cv,
    model = c("BLINK"),
    PCA.total = 0,
    file.output = T,
    group.from = 616,
    group.to=616,
    Inter.Plot=FALSE,
  )
}

