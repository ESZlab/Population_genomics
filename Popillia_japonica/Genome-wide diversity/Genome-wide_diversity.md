# POPULATION STRUCTURE AND DIVERSITY ANALYSIS USING STACKS

### <u>CITATION</u>: 
Catchen J., Hohenlohe P.A., Bassham S., Amores A., Cresko W.A. (2013). Stacks: an analysis tool set for population genomics. *Molecular Ecology*.

### <u>USEFUL LINK</u>: 
https://catchenlab.life.illinois.edu/stacks/comp/populations.php


## DESCRIPTION:
The `populations` module from Stacks was used to calculate key population genomic indices:
- Average observed heterozygosity (Ho)
- Average expected heterozygosity (He)
- Inbreeding coefficient (Fis)


## BEFORE STARTING:
A population map file is required. This file must contain:
- One individual per line (using sample names exactly as in the VCF file)
- A tab-delimited structure: `<sample_name><TAB><population_name>`

Example:

```
head -5 Pj_popmap.txt

DMR120j	South Japan
DMR122j	South Japan
DMR125j	South Japan
DMR127j	South Japan
DMR128j	South Japan
```

## SETUP: create and activate a conda environment for Stacks
```
conda create -n stack_env -c bioconda stacks && conda activate stack_env
```

## RUN STACKS `populations` MODULE

```
populations -V ~/snps_3p.vcf -O ./ -M Pj_6p_popmap.txt
```

Where:

`-V`: Input VCF file

`-O`: Output directory (current directory in this example)

`-M`: Population map file


## NOTES
- This analysis may take ~2 hours for a VCF file containing ~3,000,000 SNPs.
 - Output includes various population genetics statistics across and within populations.

# ------------------------------------------------------
# PER-POPULATION SLIDING WINDOW ANALYSIS USING VCFtools v0.1.16

### <u>CITATION</u>: 
Danecek P., Auton A., Abecasis G.R., Albers C.A., Banks E., DePristo M.A.,  Handsaker R.E., Lunter G., Marth G.T., Sherry S.T., McVean G., Durbin R. (2011).
The Variant Call Format and VCFtools. *Bioinformatics*.

## Input files:
 - VCF file: snps_3p_new.vcf
 - Window width: 5 kb (set later during sliding window step)
 - TXT files containing sample names for each population (sample names must match those in the VCF)
 - Population-specific VCF files are required for individual analyses

## STEP 1: Extract population-specific VCFs using sample name lists

#### Required format:
`vcftools --vcf <vcf_file> --keep <population_namelist.txt> --out <output_prefix> --recode`

### Commands for each population:
```
vcftools --vcf snps_3p_new.vcf --keep namelist_e_popmap/sjap.txt --out sjap --recode        # South Japan  
vcftools --vcf snps_3p_new.vcf --keep namelist_e_popmap/ncjap.txt --out ncjap --recode      # North/Centre Japan  
vcftools --vcf snps_3p_new.vcf --keep namelist_e_popmap/usca.txt --out usca --recode        # USA + Canada  
vcftools --vcf snps_3p_new.vcf --keep namelist_e_popmap/azjor.txt --out azjor --recode      # Azores: São Jorge  
vcftools --vcf snps_3p_new.vcf --keep namelist_e_popmap/azmig.txt --out azmig --recode      # Azores: São Miguel  
vcftools --vcf snps_3p_new.vcf --keep namelist_e_popmap/ittc.txt --out ittc --recode        # Italy + Ticino  
```

### Note:
 - The `--recode` flag is necessary to generate a new VCF file for each population.
- Output files will have the `.recode.vcf` extension by default (e.g., sjap.recode.vcf).
- It's recommended to rename these files to remove the `.recode` part, e.g.:
```
mv sjap.recode.vcf sjap.vcf
```

# CALCULATION OF NUCLEOTIDE DIVERSITY (π) WITH VCFtools v0.1.16

### <u>CITATION</u>: 
The Variant Call Format and VCFtools — Petr Danecek et al., *Bioinformatics*, 2011.

## DESCRIPTION:
For each population-specific VCF file, we calculate nucleotide diversity (π) using non-overlapping sliding windows of 5,000 bp, where:

```
--window-pi <int> : Computes π in windows of defined size (in bp).
                     Output file will have the suffix: ".windowed.pi"
```
## COMMANDS FOR EACH POPULATION VCF:
```
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/sjap.vcf   --window-pi 5000 --out sjap_pi_5kb      # South Japan
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/ncjap.vcf  --window-pi 5000 --out ncjap_pi_5kb     # North/Centre Japan
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/usca.vcf   --window-pi 5000 --out usca_pi_5kb      # US + Canada
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/azjor.vcf  --window-pi 5000 --out azjor_pi_5kb     # Azores: São Jorge
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/azmig.vcf  --window-pi 5000 --out azmig_pi_5kb     # Azores: São Miguel
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/ittc.vcf   --window-pi 5000 --out ittc_pi_5kb      # Italy + Ticino
```
## NOTES:
- Output files: `*.windowed.pi`
- Each line of the output includes: `CHROM, BIN_START, BIN_END, N_VARIANTS, PI`
- Suitable for downstream visualization and comparative analysis

# CALCULATION OF TAJIMA'S D (5 kb windows) WITH VCFtools v0.1.16
## DESCRIPTION:
Calculate the Tajima’s D statistic in non-overlapping windows of the specified size (in bp):

`--TajimaD <integter> `

The output file will have the suffix `*.Tajima.D`.

Input VCF files: one per population

Output: one file per population, containing Tajima’s D values in 5 kb windows
```
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/sjap.vcf   --TajimaD 5000 --out sjap_tjd_5kb     # South Japan
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/ncjap.vcf  --TajimaD 5000 --out ncjap_tjd_5kb    # North/Centre Japan
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/usca.vcf   --TajimaD 5000 --out usca_tjd_5kb     # US + Canada
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/azjor.vcf  --TajimaD 5000 --out azjor_tjd_5kb    # Azores: São Jorge
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/azmig.vcf  --TajimaD 5000 --out azmig_tjd_5kb    # Azores: São Miguel
vcftools --vcf ~/snps_pj/datasets_vcf/sep_pop_vcf/ittc.vcf   --TajimaD 5000 --out ittc_tjd_5kb     # Italy + Ticino
```
# Downstream Analysis and Visualization of Sliding Windows Analysis in R

```
## Plotting Sliding Window Analysis - Nucleotide Diversity (π)
# Load necessary libraries for data manipulation and visualization
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)
library(hrbrthemes)
library(tidyverse)

# Set working directory where result files are stored
setwd("~/SNPs_analysis/vcftools/perpopulation_SlWind/PI/results/")

# Load datasets for each population's sliding window analysis of π
sjap <- read.table("sjap_pi_5kb.windowed.pi", header=T)
ncjap <- read.table("ncjap_pi_5kb.windowed.pi", header=T)
usca <- read.table("usca_pi_5kb.windowed.pi", header=T)
azjor <- read.table("azjor_pi_5kb.windowed.pi", header=T)
azmig <- read.table("azmig_pi_5kb.windowed.pi", header=T)
ittc <- read.table("ittc_pi_5kb.windowed.pi", header=T)

# Calculate median and 5th and 95th percentiles for each population
median(sjap$PI)
quantile(sjap$PI, probs = c(0.05,0.95))
summary(sjap$PI)

# Repeat for other populations
median(ncjap$PI)
quantile(ncjap$PI, probs = c(0.05,0.95))
summary(ncjap$PI)

median(usca$PI)
quantile(usca$PI, probs = c(0.05,0.95))
summary(usca$PI)

median(azjor$PI)
quantile(azjor$PI, probs = c(0.05,0.95))
summary(azjor$PI)

median(azmig$PI)
quantile(azmig$PI, probs = c(0.05,0.95))
summary(azmig$PI)

median(ittc$PI)
quantile(ittc$PI, probs = c(0.05,0.95))
summary(ittc$PI)

# **Violin Plots for Visualization**
# Combine data from all populations into a single dataset for plotting
data <- data.frame(
  name = c(rep("South Japan", 104575), rep("North/Centre Japan", 106264), 
           rep("US and Canada", 106167), rep("Sao Jorge", 103819), 
           rep("Sao Miguel", 92596), rep("Italy and Ticino", 105736)),
  value = c(sjap$PI, ncjap$PI, usca$PI, azjor$PI, azmig$PI, ittc$PI)
)

# Create the violin plot for nucleotide diversity (π) across populations
ggplot(data, aes(x=name, y=value, fill=name)) + 
  geom_violin()

# Calculate the sample size for each population (useful for further analysis)
sample_size = data %>% group_by(name) %>% summarize(num=n()) # Optional: Used for plotting by sample size

# Enhanced violin plot with boxplot overlay, statistical significance, and custom colors
data %>%
  mutate(name = fct_relevel(name, 
                            "South Japan", "North/Centre Japan", "US and Canada", 
                            "Sao Jorge", "Sao Miguel", "Italy and Ticino")) %>%
  ggplot(aes(name, value)) +
  geom_violin(width=0.9, linewidth=0.2, aes(fill=name)) +
  geom_boxplot(width=0.1, color="grey50", alpha=0.2) +
  stat_summary(fun=median, geom='point', colour='black', size=2) +
  scale_fill_manual(values = c('#00560a','#3bff49','#000080','#ff7f00', '#ffe100','#e31a1c')) +
  scale_y_continuous() +
  stat_pvalue_manual(pairwise_comparisons, label = "p.adj.signif", y.position = "y.position", linetype = 'solid') + 
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Sliding Window Analysis - Nucleotide Diversity (π)") +
  geom_hline(yintercept=0.0, linetype="dashed", color = "black") +
  xlab('') +
  ylab("Nucleotide diversity (π)")

# **Statistical Analysis and Significance Testing**
# Load necessary libraries for statistical testing
library(ggpubr)
library(rstatix)
library(ggplot2)
library(dplyr)
library(forcats)

# Function to calculate standard error
calculate_se <- function(data) {
  n <- nrow(data)  # Number of windows
  s <- sd(data$PI) # Standard deviation of PI values
  se <- s / sqrt(n) # Standard error
  return(se)
}

# Calculate standard error for each population
se_sjap <- calculate_se(sjap)
se_ncjap <- calculate_se(ncjap)
se_usca <- calculate_se(usca)
se_azjor <- calculate_se(azjor)
se_azmig <- calculate_se(azmig)
se_ittc <- calculate_se(ittc)

# Print the standard errors for each population
cat("Standard Error for sjap:", se_sjap, "\n")
cat("Standard Error for ncjap:", se_ncjap, "\n")
cat("Standard Error for usca:", se_usca, "\n")
cat("Standard Error for azjor:", se_azjor, "\n")
cat("Standard Error for azmig:", se_azmig, "\n")
cat("Standard Error for ittc:", se_ittc, "\n")

# Perform Kruskal-Wallis test for overall significance between populations
kusk = kruskal.test(value ~ name, data = data)

# Prepare pairwise comparisons for plotting
pairwise_comparisons <- data %>%
  wilcox_test(value ~ name, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "boh", step.increase = .4)

# Set colors for each population group
mycol <- c('South Japan' = '#00560a', 'North/Centre Japan' = '#3bf149', 
           'US and Canada' = '#000080', 'Sao Jorge' = '#ff7f00', 
           'Sao Miguel' = "#ffe100", 'Italy and Ticino' = '#de1a1c')

# Ensure correct factor levels for populations
data$name <- factor(data$name, levels = c("South Japan", "North/Centre Japan", "US and Canada", 
                                          "Sao Jorge", "Sao Miguel", "Italy and Ticino"))

# Violin plot with statistical significance, color adjustments, and cleaner labeling
ggplot(data, aes(x = name, y = value)) +
  geom_violin(width = .9, linewidth = 0.2, aes(fill=name)) +
  geom_boxplot(width = 0.1, color = "grey40", alpha = 0.2) +
  stat_pvalue_manual(pairwise_comparisons, label = "p.adj.signif", y.position = "y.position", linetype = 'longdash') +
  stat_summary(fun = median, geom = "point", colour = "black", size = 2) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous() +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size = 11)) +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "black") +
  ylab("π, Nucleotide Diversity") +
  xlab("")

# **Plot Adjustment**
max_value <- max(data$value)
range_value <- diff(range(data$value))
adjustment_factor <- 0.1 * range_value  # Small adjustment factor to prevent overlap

# Increment y.position to avoid overlap in pairwise comparisons
increment <- seq(0, by = adjustment_factor / 2, length.out = nrow(pairwise_comparisons))
pairwise_comparisons$y.position <- max_value + adjustment_factor + increment

# Adjust xmin and xmax for pairwise comparisons
pairwise_comparisons <- pairwise_comparisons %>%
  mutate(
    xmin = as.numeric(factor(group1, levels = levels(data$name))) - 0,
    xmax = as.numeric(factor(group2, levels = levels(data$name))) + 0
  )

# Final plot with better formatting and saved as PDF
ggplot(data, aes(x = name, y = value)) +
  geom_violin(width = .9, linewidth = 0.2, aes(fill = name)) +
  geom_boxplot(width = 0.1, color = "grey40", alpha = 0.2) +
  stat_summary(fun = median, geom = "point", colour = "black", size = 1) +
  scale_fill_manual(values = mycol) +
  scale_y_continuous(breaks=seq(0,0.015,0.0025)) +
  stat_pvalue_manual(pairwise_comparisons, 
                     label = "p.adj.signif", 
                     y.position = pairwise_comparisons$y.position, 
                     linetype = 'solid', 
                     tip.length = 0) +
  theme_classic() +
  theme(legend.position = "none", plot.title = element_text(size = 11)) +
  ggtitle("Sliding Window Analysis - Nucleotide Diversity (π)") +
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "black") +
  ylab("π, Nucleotide Diversity") +
  xlab("")
  
  
## Plotting Sliding Window Analysis - Tajima's D

# Set working directory where result files are stored
setwd("/home/rf/Scrivania/SNPs_analysis/vcftools/perpopulation_SlWind/TAJIMA'S D")

# Load and clean datasets for each population
sjap <- read.table("sjap_tjd_5kb.Tajima.D", header=T) %>% na.omit()
ncjap <- read.table("ncjap_tjd_5kb.Tajima.D", header=T) %>% na.omit()
usca <- read.table("usca_tjd_5kb.Tajima.D", header=T) %>% na.omit()
azjor <- read.table("azjor_tjd_5kb.Tajima.D", header=T) %>% na.omit()
azmig <- read.table("azmig_tjd_5kb.Tajima.D", header=T) %>% na.omit()
ittc <- read.table("ittc_tjd_5kb.Tajima.D", header=T) %>% na.omit()

# Calculate basic statistics: median, mean, and quantiles for each population
# South Japan
median(sjap$TajimaD)
mean(sjap$TajimaD)
quantile(sjap$TajimaD, probs = c(0.05, 0.95))

# North/Centre Japan
median(ncjap$TajimaD)
mean(ncjap$TajimaD)
quantile(ncjap$TajimaD, probs = c(0.05, 0.95))

# US and Canada
median(usca$TajimaD)
mean(usca$TajimaD)
quantile(usca$TajimaD, probs = c(0.05, 0.95))

# Sao Jorge
median(azjor$TajimaD)
mean(azjor$TajimaD)
quantile(azjor$TajimaD, probs = c(0.05, 0.95))

# Sao Miguel
median(azmig$TajimaD)
mean(azmig$TajimaD)
quantile(azmig$TajimaD, probs = c(0.05, 0.95))

# Italy and Ticino
median(ittc$TajimaD)
mean(ittc$TajimaD)
quantile(ittc$TajimaD, probs = c(0.05, 0.95))

# Create dataset for plotting
data <- data.frame(
  name = c(rep("South Japan", length(sjap$TajimaD)), rep("North/Centre Japan", length(ncjap$TajimaD)), 
           rep("US and Canada", length(usca$TajimaD)), rep("Sao Jorge", length(azjor$TajimaD)), 
           rep("Sao Miguel", length(azmig$TajimaD)), rep("Italy and Ticino", length(ittc$TajimaD))),
  value = c(sjap$TajimaD, ncjap$TajimaD, usca$TajimaD, azjor$TajimaD, azmig$TajimaD, ittc$TajimaD)
)

# Plotting with ggplot2
data %>%
  mutate(name = fct_relevel(name, 
                            "South Japan", "North/Centre Japan", "US and Canada", 
                            "Sao Jorge", "Sao Miguel", "Italy and Ticino")) %>%
  ggplot(aes(name, value, fill=name)) +
  geom_violin(width=0.9, size=0.2) +  # Violin plot
  geom_boxplot(width=0.1, color="grey50", alpha=0.2) +  # Add boxplot for comparison
  stat_summary(fun=median, geom='point', colour='black', size=2) +  # Median point
  scale_fill_manual(values = c('#00560a','#3bff49','#000080','#ff7f00','#ffe100','#e31a1c')) +  # Custom colors
  scale_y_continuous(breaks=seq(-2.5, 4, 0.5)) +  # Y-axis range
  theme_classic() +  # Clean theme
  theme(
    legend.position="none",  # Hide legend
    plot.title = element_text(size=11)  # Title font size
  ) +
  ggtitle("Sliding Window Analysis of Tajima's D") +  # Plot title
  geom_hline(yintercept=0.0, linetype="dashed", color="black") +  # Reference line at y=0
  xlab('') +  # Remove x-axis label
  ylab("Tajima's D")  # Y-axis label
  
```

# LD Decay Analysis - PopLDdecay

### <u>CITATION</u>:
Zhang C, Dong SS, Xu JY, He WM, Yang TL. (2019). PopLDdecay: a fast and effective tool for linkage disequilibrium decay analysis based on variant call format files. *Bioinformatics*, 35(10):1786-1788.

## DESCRIPTION:
LD decay (R²) was estimated on subsets of 6 randomly selected individuals per population.
Subsetting was done via `.txt` namelists (see **Sliding Windows section**).

## Inspect VCF (optional check):
```
cat filename.vcf | grep '#' | tail
```

## PopLDdecay INSTALLATION:
```
git clone https://github.com/hewm2008/PopLDdecay.git 
cd PopLDdecay
chmod 755 configure && ./configure
make
mv PopLDdecay bin/
```

## RUN ANALYSIS (per population):
```
./bin/PopLDdecay -InVCF ~/snps_pj/datasets_snps/subsettato/sj.vcf -MaxDist 100 -OutStat SJ_LDdecay
perl bin/Plot_OnePop.pl -inFile SJ_LDdecay.stat.gz -output plot/SJ

./bin/PopLDdecay -InVCF ~/snps_pj/datasets_snps/subsettato/sub_ncj.vcf -MaxDist 100 -OutStat sub_NCJ_LDdecay
perl bin/Plot_OnePop.pl -inFile sub_NCJ_LDdecay.stat.gz -output plot/NCJ

./bin/PopLDdecay -InVCF ~/snps_pj/datasets_snps/subsettato/sub_usca.vcf -MaxDist 100 -OutStat sub_USCA_LDdecay
perl bin/Plot_OnePop.pl -inFile sub_USCA_LDdecay.stat.gz -output plot/USCA

./bin/PopLDdecay -InVCF ~/snps_pj/datasets_snps/subsettato/sub_azjor.vcf -MaxDist 100 -OutStat sub_AZJOR_LDdecay
perl bin/Plot_OnePop.pl -inFile sub_AZJOR_LDdecay.stat.gz -output plot/AZJOR

./bin/PopLDdecay -InVCF ~/snps_pj/datasets_snps/subsettato/sub_azmig.vcf -MaxDist 100 -OutStat sub_AZMIG_LDdecay
perl bin/Plot_OnePop.pl -inFile sub_AZMIG_LDdecay.stat.gz -output plot/AZMIG

./bin/PopLDdecay -InVCF ~/snps_pj/datasets_snps/subsettato/sub_ittc.vcf -MaxDist 100 -OutStat sub_ITTC_LDdecay
perl bin/Plot_OnePop.pl -inFile sub_ITTC_LDdecay.stat.gz -output plot/ITTC
```

## COMBINE POPULATIONS:
Prepare a `multi.list` file with two columns: 

`<path_to_stat.gz><TAB><population_code>`

Example:
```
/path/to/SJ_LDdecay.stat.gz SJ
/path/to/sub_NCJ_LDdecay.stat.gz NCJ
/path/to/sub_USCA_LDdecay.stat.gz USCA
/path/to/sub_AZJOR_LDdecay.stat.gz AZJ
/path/to/sub_AZMIG_LDdecay.stat.gz AZM
/path/to/sub_ITTC_LDdecay.stat.gz ITTC
```
```
perl bin/Plot_MultiPop.pl -inList multi.list -output ALL
```
## Edit color palette (Optional):

```
nano bin/Plot_MultiPop.pl
```
Modify colors, then re-run the above line.
