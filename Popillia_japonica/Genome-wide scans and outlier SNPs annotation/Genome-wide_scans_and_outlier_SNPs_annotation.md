# GENOME-WIDE SCANS FOR DETECTION OF OUTLIER SNPs USING PCAdapt IN R

### <u>CITATION</u>: 
Privé F, Luu K, Vilhjálmsson BJ, Blum MG (2020). Performing highly efficient genome scans for local adaptation with R package pcadapt version 4. *Molecular Biology and Evolution*, 37(7), 2153–2154.

### <u>USEFUL LINK</u>: 
https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

### NOTE:
Workflow applied to multiple population contrasts of *Popillia japonica*:
Example shown here: North/Central Japan vs USA + Canada

Other contrasts analyzed using the same pipeline include:
- USA + Canada vs São Miguel
- USA + Canada vs São Jorge
- São Miguel vs São Jorge
- USA + Canada vs Italy + Ticino

## Step 1: Convert VCF to PLINK format
```
plink --vcf vcf_filename.vcf --double-id --allow-extra-chr \
        --set-missing-var-ids @:# --make-bed --out outputname
```
## Step 2: Extract chromosome and position info from VCF using bcftools
```
bcftools query -f '%CHROM %POS\n' vcf_filename.vcf > scaffoldlist_name.txt
```
## Step 3: Run PCAdapt in R
```
library(pcadapt)
library(qvalue)
library(qqman)
library(dplyr)

# Load the PLINK file
pj_pcadapt_nu <- read.pcadapt("PCadapt/ncj_usca/ncj_usca.bed", type = "bed")

# Screeplot for PC selection
x_nu <- pcadapt(input = pj_pcadapt_nu, K = 20)
plot(x_nu, option = "screeplot")          # Full scree plot
plot(x_nu, option = "screeplot", K = 10)  # Zoomed in (10 PCs)

# Based on Cattell’s rule, select K = 2 (Patterson et al. 2006)
x_nu <- pcadapt(pj_pcadapt_nu, K = 2)

# Summary and diagnostic plots
summary(x_nu)
plot(x_nu, option = "manhattan")
plot(x_nu, option = "qqplot")
plot(x_nu, option = "scores")
hist(x_nu$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x_nu, option = "stat.distribution")

# -----------------------------
# Step 4: Identify outliers using q-values (FDR control)
# -----------------------------
alpha <- 0.001  # 0.1% FDR threshold
qval_nu <- qvalue(x_nu$pvalues)$qvalues
outliers_qval_nu <- which(qval_nu < alpha)
length(outliers_qval_nu)
print(outliers_qval_nu)

# Add chromosome and position info
chrm <- read.delim("scaffoldlist_name.txt", sep = "\t", header = FALSE)
qval_table_nu <- cbind(as.data.frame(qval_nu), chrm)
outliers_qval_nu_df <- filter(qval_table_nu, qval_nu < alpha) %>%
  rename(LocusName = V2)
head(outliers_qval_nu_df)

# Save q-value based outliers
write.table(outliers_qval_nu_df, "ncj_usca_outliers_qval.txt", sep = "\t", row.names = FALSE)

# -----------------------------
# Step 5: Identify outliers using Bonferroni correction  # chosen procedure
# -----------------------------
padj_bonf_nu <- p.adjust(x_nu$pvalues, method = "bonferroni")
outliers_bonf_nu <- which(padj_bonf_nu < alpha)
length(outliers_bonf_nu)
print(outliers_bonf_nu)

# Annotate Bonferroni-adjusted p-values
padj_bonf_table_nu <- cbind(as.data.frame(padj_bonf_nu), chrm) %>%
  rename(Scaffold = V1, LocusName = V2)
outliers_padj_bonf_nu <- filter(padj_bonf_table_nu, padj_bonf_nu < alpha)

# Save Bonferroni-based outliers
write.table(outliers_padj_bonf_nu, "ncj_usca_outliers_padj_bonf.txt", sep = "\t", row.names = FALSE)

# Clean file for use in overlaps
outliers_padj_bonf_clean_nu <- outliers_padj_bonf_nu[, c("Scaffold", "LocusName")]
write.table(outliers_padj_bonf_clean_nu,
            file = "ncj_usca_outliers_pcadapt.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# The following step has to be run after F<sub>ST</sub> analysis
# -----------------------------
# Step 6: Manhattan plot highlighting shared F<sub>ST</sub>/PCAdapt outliers
# -----------------------------
# Prepare data for qqman
padj_bonf_table_nu <- na.omit(padj_bonf_table_nu)
colnames(padj_bonf_table_nu) <- c("P", "CHR", "BP")
padj_bonf_table_nu$CHR <- as.numeric(gsub("[^0-9]", "", padj_bonf_table_nu$CHR))
padj_bonf_table_nu$SNP <- paste0(padj_bonf_table_nu$CHR, ":", padj_bonf_table_nu$BP)

# Highlight SNPs previously found in F<sub>ST</sub> sliding window analysis
highlight_snps_nu <- c("44:574340", "47:1388008", "64:1187582", 
                       "71:1412535", "159:567053", "183:75772", 
                       "207:872553", "226:359789", "255:685705", 
                       "406:350909", "445:99440", "450:76955", 
                       "463:365496", "706:9597", "815:137251", 
                       "962:6645", "1282:12067")

# Verify which highlighted SNPs are in the dataset
matching_snps_nu <- highlight_snps_nu[highlight_snps_nu %in% padj_bonf_table_nu$SNP]
print(matching_snps_nu)

# Plot Manhattan plot
manhattan(padj_bonf_table_nu, 
          chr = "CHR", 
          bp = "BP", 
          p = "P", 
          snp = "SNP",
          genomewideline = -log10(alpha), 
          suggestiveline = FALSE,
          cex.lab = 1.5,
          logp = TRUE,
          highlight = highlight_snps_nu,
          col = c("grey20"), 
          ylim = c(0, 50),
          xaxt = "n",
          xlab = "SNPs",
          main = "North/Central Japan vs USA and Canada")

abline(h = -log10(alpha), col = "red", lwd = 2.5)
```

# GENOME-WIDE SCANS FOR DETECTION OF OUTLIER SNPs with F<sub>ST</sub> SLIDING WINDOWS using VCFtools 


Example for *Popillia japonica* with: NORTH/CENTRE JAPAN vs USA + CANADA

subset the VCF file in order to have a VCF containg only individuals from contrasted populations using the appropriate namelist `.txt` file  
```
vcftools --vcf snps_3p_new.vcf --keep ncjap_usca.txt --out ncjap_usca --recode
```

## Run F<sub>ST</sub> sliding windows analysis 
Prepare population maps `.txt` files for each contrasted population 

```
vcftools --vcf ~/snps_pj/datasets_vcf/nat_pop_vcf/ncjap_usca.vcf --weir-fst-pop ~/snps_pj/datasets_vcf/namelist_e_popmap/ncjap_popmap.txt --weir-fst-pop ~/snps_pj/datasets_vcf/namelist_e_popmap/perpop_popmaps/usca_popmap.txt --fst-window-size 5000 --out _usca_fst_5kb
```

### NOTES: 
`--weir-fst-pop`--> This option is used to calculate an Weir and Cockerham’s F<sub>ST</sub> estimate. The provided file must contain a list of individuals (one individual per line) from the VCF file that correspond to one population. By default, calculations are done on a per-site basis. The output file has the suffix `*.weir.fst`.

`--fst-window-size <integer>`  --> These options can be used with `--weir-fst-pop` to do the F<sub>ST</sub> calculations on a windowed basis instead of a per-site basis. The argument specify the desired window size

## In R:

```
### Example with NCJ_VS_USCA ###
#set working directory
setwd('/home/rf/Scrivania/Revision_part2/revision_2/SW/')

#NCJ vs USCA
#load dataset - letters in scaffold names have to be removed 
ncj_us_fst <- read.table("noF/ncj_usca_fst_5kb.windowed.weir.fst", header=T)
head(ncj_us_fst)
summary(object = ncj_us_fst$MEAN_FST) 
summary(object = ncj_us_fst$WEIGHTED_FST)
#Min.      1st Qu.  Median   Mean     3rd Qu.  Max. 
#-0.073230 -0.003133  0.036539  0.072907  0.117447  0.818834 

#outliers preview
#slice_max selects rows with highest values of a variable. I this case we select the highest 0.1%
ou_ncjus <- slice_max(ncj_us_fst, prop = 0.001, WEIGHTED_FST)
write.csv(ou_ncjus, file = "ou_ncjus.csv", row.names = FALSE)

#finding outliers
ncjus_wfst <- c(ncj_us_fst$WEIGHTED_FST)
quantile(ncjus_wfst, 0.999) ##use this number as threshold --> 0.6852329 
manhattan(ncj_us_fst, 
          chr="CHROM",bp="BIN_START",
          p="WEIGHTED_FST",snp="N_VARIANTS",
          col='grey40',
          logp=FALSE,
          xlab='SNPs', ylab='Fst',
          genomewideline= 0.6852329,
          title('North Central Japan vs US + Canada'))

###### produce a .bed file
ncj_usca_fst_1 <- read.table("SW/ncj_usca_fst_5kb.windowed.weir.fst", header=T)
ncj_usca_wfst_outlier_window = ncj_usca_fst_1[which(ncj_usca_fst_1[,5]>0.6852329),]  ##get windows that have outlier Fst
ncj_usca_wfst_outlier_window_bed = subset(ncj_usca_wfst_outlier_window, select = CHROM:BIN_END) #extract bed file of 99th percentile windows
write.table(ncj_usca_wfst_outlier_window_bed, "outliers_fst/ncj_usca_wfst_outlier_window.bed", sep = "\t", col.names = F, row.names = F, quote = FALSE)
```

Then extract the list of SNPs to compare with PCadapt results (Do the same after each .bed)

## Running VCFtools to extract outlier snps in these windows:
```
vcftools --vcf ncj+usca.vcf --bed ncj_usca_wfst_outlier_window.bed --recode --out outliers_ncj_usca
```
get the list of outliers snps:
```
cat outliers_ncj_usca.vcf | grep -v "#" | cut -f1-2 > ncj_usca_outlier_snps.txt
```
Highlight those in common with PCAadapt

## Load the list of SNPs from a text file, in R:
```
snp_list_ncjus <- read.table("outliers_in_common/noF/outliers_common_ncj_usca.txt", header = FALSE, stringsAsFactors = FALSE)

# Assign column names corresponding to your data (assuming first column is CHROM and second is POSITION)
colnames(snp_list_ncjus) <- c("CHR", "POS")

# Ensure CHR columns in both datasets are in the same format
ncj_us_fst$CHROM <- as.numeric(ncj_us_fst$CHROM)  # Ensure CHR column is character type
snp_list_ncjus$CHR <- as.numeric(snp_list_ncjus$CHR)

# Initialize a logical vector to store the matches
matched_indices_ncjus <- logical(nrow(ncj_us_fst))

# Subset the F<sub>ST</sub> data to get only the matched windows
matched_windows_ncjus <- ncj_us_fst[matched_indices_ncjus, ]

# Ensure the highlight column is logical and contains TRUE/FALSE values
ncj_us_fst$highlight <- FALSE

# Loop through SNP list and update the highlight column
for (i in 1:nrow(snp_list_ncjus)) {
  matched_indices_ncjus <- which(ncj_us_fst$CHROM == snp_list_ncjus$CHR[i] &
                                  snp_list_ncjus$POS[i] >= ncj_us_fst$BIN_START &
                                  snp_list_ncjus$POS[i] <= ncj_us_fst$BIN_END)
  if (length(matched_indices_ncjus) > 0) {
    ncj_us_fst$highlight[matched_indices_ncjus] <- TRUE
  }
}

# Create a column for SNP identifiers to be used for highlighting
ncj_us_fst$snp_id <- ifelse(ncj_us_fst$highlight, paste0("SNP", seq_len(nrow(ncj_us_fst))), NA)

# Extract highlighted SNPs
highlighted_snp_ids_ncjus <- na.omit(ncj_us_fst$snp_id)

# Plot the Manhattan plot with proper highlighting
par(mar = c(5, 5, 3, 2))  # Adjusting the margins (bottom, left, top, right)
manhattan(ncj_us_fst, 
          chr = "CHROM", 
          bp = "BIN_START", 
          p = "WEIGHTED_F<sub>ST</sub>", 
          snp = "snp_id",  # Column for SNP identifiers
          highlight = highlighted_snp_ids_ncjus,  # SNPs to highlight
          col = c('grey20'),  # Grey for non-highlighted, red for highlighted
          logp = FALSE,
          xaxt = "n",
          ylim = c(0, 1.1),  # Adjust y-axis limits
          genomewideline = 0.6852329,
          xlab = "SNPs",
          ylab = expression(Weighted~F[ST]),
          cex.lab = 1.5,
          main = "North/Central Japan vs US and Canada")

abline(h = 0.6852329, col = "red", lwd = 2.5) 

```
# GENOME-WIDE SCANS FOR DETECTION OF OUTLIER SNPs USING F<sub>ST</sub> SLIDING WINDOWS with VCFtools

Example for *Popillia japonica* with: NORTH/CENTRE JAPAN vs USA + CANADA

## STEP 1: SUBSET VCF FILE FOR CONTRASTED POPULATIONS

Subset the VCF file to include only individuals from contrasted populations
```
vcftools --vcf snps_3p_new.vcf --keep ncjap_usca.txt --out ncjap_usca --recode
```


## STEP 2: RUN F<sub>ST</sub> SLIDING WINDOW ANALYSIS

Prepare two population maps `.txt` with individual IDs per population
```
vcftools --vcf ~/snps_pj/datasets_vcf/nat_pop_vcf/ncjap_usca.vcf \
         --weir-fst-pop ~/snps_pj/datasets_vcf/namelist_e_popmap/ncjap_popmap.txt \
         --weir-fst-pop ~/snps_pj/datasets_vcf/namelist_e_popmap/perpop_popmaps/usca_popmap.txt \
         --fst-window-size 5000 --out _usca_fst_5kb
```
### Notes:
`--weir-fst-pop`: Weir and Cockerham’s Fst estimate per population

`--fst-window-size`: Performs Fst calculation over windows (here, 5kb)

## STEP 3: OUTLIER DETECTION IN R

```
# Load the dataset (remove any scaffold suffixes if present)
ncj_us_fst <- read.table("noF/ncj_usca_fst_5kb.windowed.weir.fst", header = TRUE)

# Summarize FST values
summary(ncj_us_fst$MEAN_FST)
summary(ncj_us_fst$WEIGHTED_FST)

## DETECT OUTLIERS
# Preview top 0.1% of WEIGHTED_FST values
ou_ncjus <- slice_max(ncj_us_fst, prop = 0.001, WEIGHTED_FST)
write.csv(ou_ncjus, file = "ou_ncjus.csv", row.names = FALSE)

# Define threshold based on 99.9th percentile
ncjus_wfst <- ncj_us_fst$WEIGHTED_FST
quantile(ncjus_wfst, 0.999)  # 0.6852329 keep this number for further use

# Manhattan plot
manhattan(ncj_us_fst, 
          chr = "CHROM", bp = "BIN_START",
          p = "WEIGHTED_FST", snp = "N_VARIANTS",
          col = "grey40", logp = FALSE,
          xlab = "SNPs", ylab = "Fst",
          genomewideline = 0.6852329,
          main = "North/Central Japan vs USA + Canada")

# EXPORT OUTLIER WINDOWS TO BED FILE
ncj_usca_fst_1 <- read.table("SW/ncj_usca_fst_5kb.windowed.weir.fst", header = TRUE)
ncj_usca_wfst_outlier_window <- ncj_usca_fst_1[ncj_usca_fst_1[,5] > 0.6852329, ]
ncj_usca_wfst_outlier_window_bed <- subset(ncj_usca_wfst_outlier_window, select = CHROM:BIN_END)
write.table(ncj_usca_wfst_outlier_window_bed, "outliers_fst/ncj_usca_wfst_outlier_window.bed", 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
```
### In terminal:
```
# Extract outlier SNPs using VCFtools
vcftools --vcf ncj+usca.vcf --bed ncj_usca_wfst_outlier_window.bed --recode --out outliers_ncj_usca

# Create a list of outlier SNP positions
cat outliers_ncj_usca.vcf | grep -v "#" | cut -f1-2 > ncj_usca_outlier_snps.txt
```
## STEP 4: COMPARE WITH PCAdapt OUTLIERS in R
```
# Load SNPs in common with PCAdapt outliers
snp_list_ncjus <- read.table("outliers_in_common/noF/outliers_common_ncj_usca.txt", 
                             header = FALSE, stringsAsFactors = FALSE)
colnames(snp_list_ncjus) <- c("CHR", "POS")

# Ensure CHROM columns are numeric
ncj_us_fst$CHROM <- as.numeric(ncj_us_fst$CHROM)
snp_list_ncjus$CHR <- as.numeric(snp_list_ncjus$CHR)

# Initialize and assign highlight column
ncj_us_fst$highlight <- FALSE

# Mark overlapping windows
for (i in 1:nrow(snp_list_ncjus)) {
  match_idx <- which(ncj_us_fst$CHROM == snp_list_ncjus$CHR[i] &
                     snp_list_ncjus$POS[i] >= ncj_us_fst$BIN_START &
                     snp_list_ncjus$POS[i] <= ncj_us_fst$BIN_END)
  if (length(match_idx) > 0) {
    ncj_us_fst$highlight[match_idx] <- TRUE
  }
}

# Assign SNP IDs to highlight
ncj_us_fst$snp_id <- ifelse(ncj_us_fst$highlight, paste0("SNP", seq_len(nrow(ncj_us_fst))), NA)
highlighted_snp_ids_ncjus <- na.omit(ncj_us_fst$snp_id)

# PLOT MANHATTAN WITH HIGHLIGHTED SNPs

par(mar = c(5, 5, 3, 2))  # Set plot margins
manhattan(ncj_us_fst, 
          chr = "CHROM", 
          bp = "BIN_START", 
          p = "WEIGHTED_FST", 
          snp = "snp_id", 
          highlight = highlighted_snp_ids_ncjus,
          col = "grey20", 
          logp = FALSE,
          xaxt = "n",
          ylim = c(0, 1.1),
          genomewideline = 0.6852329,
          xlab = "SNPs",
          ylab = expression(Weighted~F[ST]),
          cex.lab = 1.5,
          main = "North/Central Japan vs USA + Canada")

abline(h = 0.6852329, col = "red", lwd = 2.5)
```
# GENOME-WIDE OUTLIER SNPs: INTERSECTIONS BETWEEN F<sub>ST</sub> AND PCADAPT RESULTS VENN DIAGRAM ANALYSIS USING R
This step allows to identify SNPs consistently identified as outliers by both F<sub>ST</sub>-based and PCAdapt-based approaches.

The procedure was replicated for all population contrasts.
```
# Required packages
library(gplots)
library(ggvenn)


# Load SNPs identified as outliers by F<sub>ST</sub> analysis
ncjap_usca_fst <- read.table("outliers_fst/venn/ncj_usca_outlier_snps.txt", header = FALSE)
colnames(ncjap_usca_fst) <- "SNP"

# Load SNPs identified as outliers by PCAdapt
ncjap_usca_pcadapt <- read.table("outliers_pcadapt/venn/ncj_usca_outliers_pcadapt.txt", header = FALSE)
colnames(ncjap_usca_pcadapt) <- "SNP"

# Build list of SNPs from both analyses
NCJ_USCA <- list(
  FST = unique(ncjap_usca_fst$SNP),
  PCAdapt = unique(ncjap_usca_pcadapt$SNP)
)

# Generate Venn diagram using gplots package
venn_output <- venn(NCJ_USCA)
print(venn_output)

# Extract and save the list of SNPs identified as outliers by both methods
ncj_usca_outliers_common <- attr(venn_output, "intersections")$`FST:PCAdapt`
write.table(
  ncj_usca_outliers_common,
  file = "outliers_in_common/outliers_common_ncj_usca.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Generate and visualize Venn diagram using ggvenn for improved aesthetics
ggvenn(NCJ_USCA, fill_color = c("#0073C2FF", "#EFC000FF")) +
  ggtitle("North/Central Japan vs US + Canada")

# Note: This procedure was repeated for each population pair considered in the genome-wide outlier SNP scan.

```

# SNPs ANNOTATION WITH SnpEff v.5.1 and SnpSift v.5.1

### <u>CITATION</u>: 
Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.". *Fly*, 6(2): 80-92.

Cingolani P, Patel VM, Coon M, Nguyen T, Land SJ, Ruden DM, Lu X. (2021). Using Drosophila melanogaster as a Model for Genotoxic Chemical Mutational Studies with a New Program, SnpSift. *Front Genet*. 15;3:35.

### <u>USEFUL LINKS</u>:
https://pcingola.github.io/SnpEff/se_build_db/

https://github.com/Elahep/BMSB-popgenomics/tree/main/SNP_annotation

snpEff is a software used for SNPs annotation and prediction of their effects on genes and proteins.

snpEff needs a database to perform genomic annotations. There are pre-built databases for thousands of genomes.

However, since we study a non-model species, the annotation information of our studied species is not available in snpEff, and we need to manually add our database to the program.

## Building SnpEff database for non-model species - *Popillia japonica*
We need to customize the config file with genomic information of our species of interest. 

It is recommended to perform the analysis in the folder where snpEff is installed

### Find the installation directory of snpEff and change to that directory
```
to_snpeff=$(find ~ -name snpEff)
cd $to_snpeff
```
### Create necessary directories
```
mkdir -p data/mygenome
mkdir -p data/genome
```

### Find the snpEff config file and copy it to the current directory
```
snpeff_contig=$(find ~ -name snpEff.config)
cp $snpeff_contig data/
```

### Find snpEff.jar and copy it as well
```
snpeff_jar=$(find ~ -name snpEff.jar)
cp snpeff_jar data/
```

In the `mygenome` folder, ensure you have:
- `genes.gff` file
- `cds.fa` file
- `protein.fa` file 

In the `genomes` folder, ensure you have:

- `mygenome.fa` file (the reference genome in fasta format)

### Now, check the config file using nano
```
nano data/snpEff.config
```
Make sure the path to the data directory is correctly set:

```
data.dir = /home/user/miniconda3/pkgs/snpeff-5.1-hdfd78af_2/share/snpeff-5.1-2/data
```

Add your species' genome details under the `Databases & Genomes` section

Add the genome entry for your species (e.g., *Popillia japonica*) in the following format:


```
# Databases & Genomes
#
# One entry per genome version. 
#
# For genome version 'ZZZ' the entries look like
#	ZZZ.genome              : Real name for ZZZ (e.g. 'Human')
#	ZZZ.reference           : [Optional] Comma separated list of URL to site/s where information for building ZZZ database was extracted.
#	ZZZ.chrName.codonTable  : [Optional] Define codon table used for chromosome 'chrName' (Default: 'codon.Standard')
#
#-------------------------------------------------------------------------------

# my genome
mygenome.genome : mygenome             ##here
```
Once the `.config `file is ready, we can build the snpEff database.

### BUILDING snpEff DATABASE WITH `.gff3` FILE AND `-noGenome` OPTION
```
java -jar snpEff.jar build -gff3 -noGenome -v mygenome
snpEff build -c snpEff.config -gff3 -noGenome -v mygenome > snpEff.stdout 2> snpEff.stderr

# Check for errors in snpEff stderr or stdout logs
cat snpEff.stderr
cat snpEff.stdout
```
## RUNNING SnpEff FOR SNP ANNOTATION
Now that the database is built, we can use snpEff to annotate SNPs from our VCF file.

Example for *Popillia japonica* with: NORTH/CENTRE JAPAN vs USA + CANADA.
This was done for each comparison in the study.

### Extract outlier SNPs from the main VCF file

```
vcftools --vcf ../PCadapt/ncj_usca/ncjap+us.vcf --snps ../outliers_in_common/outliers_common_ncj_usca.txt --recode --recode-INFO-all --out ncj_usca_annot
```

### Run SNP annotation using snpEff
```
java -Xmx4g -jar snpEff.jar -c snpEff.config mygenome out_NCJ_USCA.vcf > USCA_AZJ/annot_NCJ_USCA.vcf
```

After each run, three output files are produced:
- `annotated.vcf `with outlier SNPs and their annotations (contains an ANN field with the properties of each SNP)
- `summary.html` (guide for interpretation of annotations: https://pcingola.github.io/SnpEff/se_outputsummary/)
- `gene_list.txt` with the list of annotated genes

To avoid overwriting files in subsequent contrast analyses, move these files into a dedicated folder (e.g., NCJvsUSCA).

### Search for annotated proteins in the gene list (example for genes g15733.t1 and g15734.t1)
```
cat genes.gff | grep -F 'g15733.t1'
cat genes.gff | grep -F 'g15734.t1'
```

## Extract detailed annotation information for the SNPs with snpSift

### Install snpSift
```
conda create -n snpsift_env -c bioconda snpsift -y && conda activate snpsift_env
```

### Extract relevant annotation fields and generate a neat table
```
cat annot_results_gff/annot_NCJ_USCA.vcf \
  | ~/miniconda3/pkgs/snpeff-5.1-hdfd78af_2/share/snpeff-5.1-2/scripts/vcfEffOnePerLine.pl \
  | java -Xmx8g -jar ~/miniconda3/pkgs/snpsift-5.1-hdfd78af_0/share/snpsift-5.1-0/SnpSift.jar extractFields - \
    CHROM \
    POS \
    ID \
    REF \
    ALT \
    AF \
    "ANN[*].ALLELE" \
    "ANN[*].EFFECT" \
    "ANN[*].IMPACT" \
    "ANN[*].GENE" \
    "ANN[*].BIOTYPE" \
    "ANN[*].HGVS_C" \
    "ANN[*].HGVS_P" \
  > outlierSNPs_NCJ_USCA.annot.sift.txt

```
# Calculate F<sub>ST</sub> values under neutral evolution using Fastsimcoal 2  

To evaluate possible effects of non-equilibrium dynamics in shaping the observed F<sub>ST</sub> values, expected F<sub>ST</sub> values at different contrasts were obtained through a simulation under neutral evolution within the best model identified in the third step of demographic analysis (see STEP-BY-STEP DEMOGRAPHIC INFERENCE WITH FASTSIMCOAL 2).


### <u>CITATION</u>: 

Excoffier, L. and M. Foll. 2011. fastsimcoal: a continuous-time coalescent simulator of genomic diversity under arbitrarily complex evolutionary scenarios. Bioinformatics 27: 1332-1334.

Excoffier, L., Dupanloup, I., Huerta-Sánchez, E., and M. Foll (2013) Robust demographic inference from genomic and SNP data. PLOS Genetics 9(10):e1003905.

### <u>LINK to Fastsimcoal manual</u>:

Fastsimcoal manual: https://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal28.pdf

## File preparation
The new simulated dataset can include a number of SNPs tailored to the analysis. In this case, a dataset of 295,396 unlinked SNPs was simulated, matching the original dataset.
To proceed, the best maxL.par file—an output of the demographic analysis performed with fastsimcoal2—must be selected. This involves first identifying the best subset, defined as the one with the highest DeltaAIC between the top two models. Within this subset, the best model is determined, after which the individual runs are evaluated. The optimal run is the one with the highest Estimated Likelihood, i.e., the value closest to the Observed Likelihood (which remains constant across models).
Once the best run of the best model within the best subset is identified, the corresponding maxL.par file should be selected and its final lines modified as described below.

```
//Number of independent loci [chromosomes]
150000000 0                                                                      
//Per chromosome: Number of linkage blocks
1
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
DNA 100 0 2.1e-9

# Notes:
# 150000000 is the number of independent DNA sequences to simulate, which has to be much higher than the number of SNPs to simulate. It is possible to start with lower numbers and adjust it to reach the desired dataset
# 1 depends on the number of blocks below
# DNA is the type of data to simulate
# 100 is the simulated sequences length
# 0 is the recombination rate
# 2.1e-9 is the chosen mutation rate
```
                
## Run the Fastsimcoal 2 command
### In terminal

Once the file is ready, it is possible to run the following command

```
./fsc28 -i modifiedMaxLpar.par -n 1 -G -g -s 295396 -C12 -B12       
```

### Notes:
`-i` specifies the par file

`-n` specifies the number of files to simulate

`-G` generates a genotype table (extension .geno)

`-g` generates an arlequin project file (extension .arp)

`-s` specifies the number of SNPs to simulate

## File conversion
Once genotype table (`.geno`) and arlequin project file with the desired number of SNPs are produced, we can convert one of them to VCF
Note that the number of samples per population that you get in this simulated SNP file is the same as your projection numbers in your initial mSFS. 

The `.arp` file can be converted to VCF format using tools such as arp2vcf, which converts fastsimcoal2-simulated DNA data into VCF format (https://rdrr.io/github/dinmatias/reconproGS/man/arp2vcf.html), or other publicly available methods.

In the present workflow, the `.geno` file was converted to VCF format using an in-house Python script, which also assigned new names to the simulated samples. The script is available under the name geno_to_vcf.py.

## Sliding Windows F<sub>ST</sub> analysis using VCFtools

Based on the newly simulated dataset, updated namelists and population maps can be generated.

Once all necessary files have been prepared, a sliding window analysis can be performed using VCFtools, following the same approach described for the main analysis:
1) Subset the VCF file for the populations being contrasted
2) Run the F<sub>ST</sub> sliding window analysis

### In R 

Example for simulated and subsetted dataset with: NORTH/CENTRE JAPAN vs USA + CANADA

```
Load the dataset (remove any scaffold suffixes if present)
ncj_us_simulated_fst <- read.table("noF/simulatated_ncj_usca_fst_5kb.windowed.weir.fst", header = TRUE)

ncjus_simulated_wfst <- ncj_us_simulated_fst$WEIGHTED_FST
mean(ncjus_simulated_wfst)
quantile(ncjus_simulated_wfst, 0.999)
```
F<sub>ST</sub> values expected under neutral evolution can now be compared to the observed values.
