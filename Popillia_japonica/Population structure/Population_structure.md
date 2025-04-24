# PRINCIPAL COMPONENT ANALYSIS (PCA) using PLINK v1.90b6.1

### <u>CITATION</u>:
Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MA, Bender D, Maller J, Sklar P, de Bakker PI, Daly MJ, Sham PC. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. *Am J Hum Genet*, 81(3), 559–575.

## Run PCA with PLINK
```
plink --vcf snps_3p_unl.vcf \
      --double-id \
      --allow-extra-chr \
      --set-missing-var-ids @:# \
      --make-bed \
      --pca \
      --out pca_pj
```
      
## PCA PLOT (PC1 vs PC2) using R
```
library(ggplot2)

# Read PCA results
eigenval <- scan("pca_pj.eigenval")
eigenvec <- read.table("pca_pj.eigenvec", header = FALSE)

# Extract only eigenvectors
eigenvectors <- eigenvec[, 3:ncol(eigenvec)]

# Population assignment (adjust if sample size changes)
populations <- c(rep("SJ", 6), rep("NCJ", 15), rep("USCA", 21), 
                 rep("AZJ", 10), rep("AZM", 10), rep("ITTC", 21))

# Adding population information to PCA results for plotting
df <- data.frame(ID = eigenvec[, 1],
                 Population = populations,
                 PC1 = eigenvectors[, 1],
                 PC2 = eigenvectors[, 2])

# Plot
ggplot(df, aes(x = PC1, y = PC2, colour = Population)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("AZJ" = '#ff9f00',
                                "AZM" = '#ffe100',
                                "ITTC" = '#e31a1c',
                                "NCJ" = '#3bff49',
                                "SJ" = '#00560a',
                                "USCA" = "#000080")) +
  labs(x = paste0("PC1 (", round(eigenval[1] / sum(eigenval) * 100, 1), "%)"),
       y = paste0("PC2 (", round(eigenval[2] / sum(eigenval) * 100, 1), "%)")) +
  theme_classic() +
  theme(legend.title = element_blank())
```

# ADMIXTURE ANALYSIS USING LEA PACKAGE in R

### <u>CITATION</u>:
Frichot E, Mathieu F, Trouillon T, Bouchard G, François O. (2014). Fast and efficient estimation of individual ancestry coefficients. *Genetics*, 196(4): 973–83.

Frichot E, François O. (2015). LEA: An R package for landscape and ecological association studies. *Methods Ecol Evol*, 6: 925–929.

## In R:
```
library(LEA)

# Convert VCF to .geno format (using unlinked SNPs)
vcf2geno(input.file = "snps_3p_unl.vcf", output.file = "Pj.geno")

# Run snmf for K from 1 to 10, with 10 repetitions each
projectalpha <- snmf("Pj.geno", K = 1:10, repetitions = 10, entropy = TRUE, CPU = 64, project = "new")

# Plot cross-entropy criterion to choose the best K (lowest value)
pdf("cross_ent_alphadefault.pdf")   
plot(projectalpha, col = "maroon4", pch = 19, cex = 1.2)
dev.off()

# Identify best runs (lowest cross-entropy) for each K
best2 <- which.min(cross.entropy(projectalpha, K = 2))
best3 <- which.min(cross.entropy(projectalpha, K = 3))
best4 <- which.min(cross.entropy(projectalpha, K = 4))
best5 <- which.min(cross.entropy(projectalpha, K = 5))
best6 <- which.min(cross.entropy(projectalpha, K = 6))
best7 <- which.min(cross.entropy(projectalpha, K = 7))
best8 <- which.min(cross.entropy(projectalpha, K = 8))
best9 <- which.min(cross.entropy(projectalpha, K = 9))
best10 <- which.min(cross.entropy(projectalpha, K = 10))

# PLOT RESULTS USING pophelper

library(pophelper)

# Load Q files from best runs (stored in "All_Qfiles/" folder)
sfiles <- list.files(path = "All_Qfiles/", full.names = TRUE)
slist <- readQ(files = sfiles)

# Plot admixture barplots for each K (excluding K=1)
plotQ(qlist = slist[2], imgtype = "pdf", height = 1.5,
      clustercol = c("#00560a", "navy"), dpi = 1200, exportpath = "./pophelper")
plotQ(qlist = slist[3], imgtype = "pdf", height = 1.5,
      clustercol = c("#ffe100", "navy", "#00560a"), dpi = 1200, exportpath = "./pophelper")
plotQ(qlist = slist[4], imgtype = "pdf", height = 1.5,
      clustercol = c("#3bff49", "#00560a", "#ffe100", "navy"), dpi = 1200, exportpath = "./pophelper")
plotQ(qlist = slist[5], imgtype = "pdf", height = 1.5,
      clustercol = c("#00560a", "#3bff49", "#ffe100", "#e31a1c", "navy"), dpi = 1200, exportpath = "./pophelper")
plotQ(qlist = slist[6], imgtype = "pdf", height = 1.5,
      clustercol = c("#ffe100", "#3bff49", "#00560a", "#ff7f00", "#e31a1c", "navy"), dpi = 1200, exportpath = "./pophelper")
plotQ(qlist = slist[7], imgtype = "pdf", height = 1.5,
      clustercol = c("pink1", "#00560a", "navy", "#3bff49", "#e31a1c", "#ff7f00", "#ffe100"), dpi = 1200, exportpath = "./pophelper")
plotQ(qlist = slist[8], imgtype = "pdf", height = 1.5,
      clustercol = c("#ffe100", "#00560a", "lightskyblue1", "#ff7f00", "#e31a1c", "navy", "pink1", "#3bff49"), dpi = 1200, exportpath = "./pophelper")
plotQ(qlist = slist[9], imgtype = "pdf", height = 1.5,
      clustercol = c("hotpink4", "lightskyblue1", "#3bff49", "#00560a", "#ffe100", "pink1", "#e31a1c", "navy", "#ff7f00"), dpi = 1200, exportpath = "./pophelper")
plotQ(qlist = slist[2], imgtype = "pdf", height = 1.5, 
      clustercol = c("lightskyblue1", "#00560a", "#20B2AA", "navy", "#ffe100", "pink1", "hotpink4", "#e31a1c", "#ff7f00", "#3bff49"), dpi = 1200, exportpath = "./pophelper") 

```
# MAXIMUM LIKELIHOOD PHYLOGENY INFERENCE USING IQ-TREE v2.0.3

### <u>CITATION</u>:
Minh B.Q. et al. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution*, 37(5), 1530–1534.

## Conversion of VCF to PHYLIP file

### <u>CITATION</u>:
Ortiz, E.M. (2019). vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. 
### <u>USEFUL LINK</u>:
https://github.com/edgardomortiz/vcf2phylip

## Step 1: Convert SNPs in VCF format to PHYLIP format
```
python vcf2phylip.py -i snps_3p_unl.vcf.gz --output-folder .
```
## Step 2: Infer phylogeny using IQ-TREE
```
iqtree -s snps_3p_unl.min4.phy -m MFP -B 2000 --seqtype DNA
```

### Notes:
- The `-m MFP` option performs ModelFinder to select the best-fit model automatically.
- The `-B 2000` option enables ultrafast bootstrap with 2000 replicates.
- The best-fit model will be reported in the `.log` file.


# CALCULATION OF PAIRWISE Fst USING StAMPP PACKAGE IN R

### <u>CITATION</u>:
Pembleton, L., Cogan, N., & Forster, J. (2013). StAMPP: an R package for calculation of genetic differentiation and structure of mixed-ploidy level populations. *Molecular Ecology Resources*, 13, 946–952.

## In R:

```
# Load required packages
library(vcfR)        # For reading and handling VCF files
library(adegenet)    # For converting to genlight objects
library(StAMPP)      # For calculating pairwise Fst values
library(ggplot2)     # (Optional) For heatmap visualization
library(corrplot)    # (Optional) For correlation-style plots
library(reshape2)    # (Optional) For melting matrices into long-format data frames

###############################
## STEP 1: Prepare input data

# Read VCF file containing SNP data
snp_vcf2 <- read.vcfR("snps_3p_unl.vcf")

# Read population map file
# NOTE: The population map should be a two-column text file without header,
# where column 1 is the sample ID and column 2 is the population name.
pop.data2 <- read.table("Pj_6p_popmap.txt", header = FALSE)

# Convert VCF to genlight object (used by adegenet and StAMPP)
gl.snp2 <- vcfR2genlight(snp_vcf2)

# Assign population information to the genlight object
pop(gl.snp2) <- pop.data2$V2

# (Optional) Perform a PCA to check population structure visually
snp.pca2 <- glPca(gl.snp2, nf = 10)

###############################
## STEP 2: Calculate pairwise Fst values between populations

# Perform pairwise Fst calculation using StAMPP
# This function also returns p-values and confidence intervals via bootstrapping
Pjap_Fst <- stamppFst(gl.snp2, nboots = 100, percent = 95, nclusters = 6)

# Extract Fst values and p-values from the output
Fst <- Pjap_Fst$Fsts
pFst <- Pjap_Fst$Pvalues

# Save results to text files
write.table(Fst, "PJ_Fst.txt", sep = "\t", quote = FALSE)
write.table(pFst, "PJ_Fst_pvalue.txt", sep = "\t", quote = FALSE)

###############################
## STEP 3: Visualize Fst matrix with a heatmap

# Read the saved Fst matrix as a numeric matrix
PJ_Fst <- as.matrix(read.table("PJ_Fst.txt", header = TRUE, row.names = 1))

# Convert matrix into long format for ggplot2
PjapFs <- melt(PJ_Fst, na.rm = TRUE)

# Check summary statistics to determine the appropriate color scale
summary(PjapFs$value)

# Generate heatmap (optional)
ggplot(data = PjapFs, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#ffd60a", mid = "#4e9de6", high = "#001d3d",       # Color palette
    midpoint = 0.13, limit = c(0.030, 0.220), space = "Lab",  # Adjust according to your data range
    name = "Pairwise Fst"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  coord_fixed()
```

# SCALED COVARIANCE MATRIX OF POPULATION ALLELE FREQUENCIES (Ω) USING BayPass v2.3 - Core Model

### <u>CITATION</u>:
Gautier M (2015). Genome-Wide Scan for Adaptive Differentiation and Association Analysis with population-specific covariables. *Genetics*, 201(4):1555–1579.

### <u>USEFUL LINKS</u>:
https://www1.montpellier.inra.fr/CBGP/software/baypass/


### Download `reshaper_baypass.py`:
https://gitlab.com/YDorant/Toolbox/-/blob/master/reshaper_baypass.py

## STEP 1: Prepare input files
Convert unlinked SNP VCF file to the BayPass .geno format using the `reshaper_baypass.py` script by Yann Dorant.

The population map file should be a two-column text file with `sample IDs` and `population labels`

## Command:
```
python ./reshaper_baypass.py ../snps_3p_unl_new.vcf ../pj_popmap.txt baypass_pj.geno
```

#### DESCRIPTION OF THE `.geno` FILE:
- Rows: SNPs
- Columns: 2 × number of populations (2 per population: one for ref allele count, one for alt allele count)
- Within each population's pair of columns:
  * First = number of reference allele copies
  * Second = number of alternative allele copies
  * "0" = no copies (e.g., missing data)

## STEP 2: Run the BayPass core model
This will compute the scaled covariance matrix of allele frequencies (Ω)
```
baypass -gfile baypass_pj.geno -outprefix pj
```
This will produce several output files, including:
- `pj_mat_omega.out` containing scaled covariance matrix (Ω)

## STEP 3: Visualize the Ω matrix as a correlation heatmap in R

```
# Load required R libraries
library(corrplot)
library(RColorBrewer)

# Read the matrix
omega <- as.matrix(read.table("pj_mat_omega.out"))

# Assign meaningful row and column names (order must match your populations)
colnames(omega) <- c('SJ', 'NCJ', 'USCA', 'AZJOR', 'AZMIG', 'ITTC')
rownames(omega) <- c('SJ', 'NCJ', 'USCA', 'AZJOR', 'AZMIG', 'ITTC')

# Convert the covariance matrix to a correlation matrix
cor.mat <- cov2cor(omega)

# Plot correlation matrix using corrplot
corrplot(
  cor.mat,
  method = "color",
  col = colorRampPalette(brewer.pal(11, "RdBu"))(200),  # Red-Blue diverging palette
  type = "upper",
  order = "hclust",     # Clustering order can help interpret patterns
  addCoef.col = "black",# Add correlation coefficients
  tl.col = "black",     # Text label color
  tl.srt = 45           # Rotate axis labels
)
```

