# STEP-BY-STEP DEMOGRAPHIC INFERENCE WITH FASTSIMCOAL 2

### <u>CITATION</u>: 
Excoffier, L. and M. Foll. 2011. fastsimcoal: a continuous-time coalescent simulator of genomic diversity under arbitrarily complex evolutionary scenarios. *Bioinformatics* 27: 1332-1334.

Excoffier, L., Dupanloup, I., Huerta-Sánchez, E., and M. Foll (2013) Robust demographic inference from genomic and SNP data. *PLOS Genetics* 9(10):e1003905.

### <u>USEFUL LINK</u>: 
https://github.com/isaacovercast/easySFS. 

## DESCRIPTION: 
This analysis follows a step-by-step demographic inference procedure.

Populations are added at each step to infer their origin based on SNP subsetting.

## STEP 1: Prepare Input Files
Start by extracting header and SNP data from VCF files.
```
grep "^#" *.vcf > temp_header.txt
grep -v "^#" *.vcf > temp_snps.txt
```
Create 5 independent random subsets of 30'000 SNPs each (choose the number of SNPs based on your data)
```
for i in {1..5}; do
  echo $i
  shuf -n 30000 temp_snps.txt | cat temp_header.txt - > subset_all_${i}.vcf
done   
```
## STEP 2: Subsetting Populations
For each subset, use vcftools to remove populations not included in the current step of the demographic model.
The example uses the populations for step 2.
```
for i in {1..5}; do
  vcftools --vcf subset_all_${i}.vcf \
           --keep namelist_step2.txt \
           --out subset_${i}_step2 \
           --recode
done
```
## STEP 3: VCF to SFS Conversion Using `easySFS.py`
Convert VCF file to SFS format. First, preview the data to check projections for each population.
```
easySFS.py -i snps_3p_unl_new.vcf -p new_popmap.txt --preview  # Preview with the full VCF and population map
```
The preview file shows the number of samples per projection and segregating sites, which helps in balancing between segregating sites and sample size. 

Once projections are reviewed, specify the projections for each population and generate SFS.

```
easySFS.py -i snps_3p_unl_new.vcf -p new_popmap.txt --proj 12,12,18,18,12  ## Projection values for each population
```
The command generates the SFS file and creates a folder named 'fastsimcoal2' with the resulting files.

The multidimensional SFS file in our case is named: `snps_3p_unl_new_MSFS.obs`

Do the same with all subsetted VCF files

## STEP 4: Fastsimcoal2 Simulation Setup
Prepare directories for each model to be tested.
`.tpl` and `.est` files have to be prepared according to Fastsimcoal manual: https://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal28.pdf.

In each directory, we need `.tpl`, `.est`, and `MSFS.obs` files. Ensure the names match the model directory.

Dicrectories and files must be arranged with the main folder corresponding to the source of introduction (e.g. Sao Miguel - AZM) where there are stored all the subset folders (e.g. I-V) each containing all the projected models (e.g. modelA_1_1 to modelA_3_1). Each model folder must contain the `.tpl`, `.obs`, `.est` , `modelA_.sh` file and a `loop.sh` script as the following example: 
```
1.AZM_Source/
├── subset_I
├── subset_II
├── subset_III
├── subset_IV
├── subset_V
   ├── modelA_1_1
   ├── modelA_2_1
   ├── modelA_3_1
      └── modelA_1_1.est
      └── modelA_1_1.tpl
      └── modelA_1_1_MSFS.obs
      └── model_.sh
      └── loop.sh
		   
```
where `model_.sh` is the main fastsimcoal2 command:
```
fsc27093 -t *.tpl -e *.est -m -n500000 -L50 -c12 -B12 --removeZeroSFS -M --multiSFS
```
and `loop.sh` will perform the loop of the 50 runs:
```
for j in {1..50}; do 
  mkdir run$j # create the runX folder
  cp *_* run$j # copy all the tpl, obs, est and sh file within the newly created run folder
  cd run$j # move in the folder
  bash ./model_.sh # launch the fastsimcoal2 script
  cd .. # move back and re-start the loop
done
```

### SUGGESTION: 
Quickly run the code with `-L 1` and `-n10000` to generate a `.par` file for model verification.

Then, convert the `.par` file to a `.pdf` using the script available at: https://cmpg.unibe.ch/software/fastsimcoal26/parFileInterpreter.html.

## STEP 5: Model Selection and Parameter Estimation
After running Fastsimcoal2, we analyze the results and select the best model according to model comparison with the Akaike Information Criterion (AIC). 

In R:
```
# Load necessary libraries
library(ggplot2)    # For plotting
library(gridExtra)  # For arranging multiple plots
library(qpcR)       # For AIC weights calculation (akaike.weights)
library(ggthemes)   # Additional themes for ggplot2 (not explicitly used here but often useful)
library(dplyr)      # For data manipulation

# ---- Define helper function ----
# Function to calculate AIC differences (deltaAIC), relative likelihoods, and Akaike weights
process_group <- function(df) {
  w <- akaike.weights(x = df$AIC)  # Compute Akaike weights using qpcR
  
  # Create dataframes for deltaAIC, relative likelihoods, and Akaike weights
  deltaAIC_results <- data.frame(modello = df$modello, deltaAIC = w$deltaAIC)
  relLL_results <- data.frame(modello = df$modello, relative_Lhoods = w$rel.LL)
  Akaikeweights_results <- data.frame(modello = df$modello, Akaike_weights = w$weights)
  
  # Merge these new columns into the original dataframe
  df <- merge(df, deltaAIC_results, by = 'modello')
  df <- merge(df, relLL_results, by = 'modello')
  df <- merge(df, Akaikeweights_results, by = 'modello')
  
  # Sort models by AIC (ascending order)
  df <- df[order(df$AIC),]
  
  return(df)
}

# ---- Set working directory and read input files ----
setwd('/home/user/AZM_Source/')

# Find all files ending with 'bestlhoods' recursively
samples <- list.files(path = ".", pattern = "bestlhoods$", recursive = TRUE, full.names = TRUE)

# Initialize an empty list to collect dataframes
tables <- list()

# Loop through each file
for (i in 1:length(samples)) {
  model <- read.table(samples[i], header = TRUE)  # Read table
  
  # Extract 'subset' info (e.g., subset_I, subset_II, etc.)
  subset_match <- regmatches(samples[i], gregexec('subset_[I-V]*', samples[i]))
  model$subset <- paste(unlist(subset_match)[1], collapse = " ")  # Add as a new column
  
  # Extract 'run' number (e.g., run1, run2, etc.)
  run_match <- regmatches(samples[i], gregexec('run[0-9]*', samples[i]))
  model$run <- paste(unlist(run_match)[1], collapse = " ")  # Add as a new column
  
  # Extract 'model' name (e.g., modelA.1, modelB.2, etc.)
  model_match <- regmatches(samples[i], gregexec('model[A-Z].[0-9]*', samples[i]))[1]
  model$modello <- paste(unlist(model_match)[1], collapse = " ")  # Add as a new column
  
  # Count the number of model parameters by matching specific column names
  parameters <- sum(
    grepl("^N_BOT", colnames(model)),
    grepl("^Admix", colnames(model)),
    grepl("^ENDBOT", colnames(model)),
    grepl("^DURBOT", colnames(model)),
    grepl("^N_POP", colnames(model))
  )
  
  model$k <- parameters  # Add number of parameters as a new column
  
  # Keep only important columns
  model <- model[, c("MaxEstLhood", "MaxObsLhood", "subset" , "modello", 'run', 'k')]
  
  # Calculate AIC (corrected for log base)
  model$AIC <- min(2 * parameters - 2 * (model$MaxEstLhood / log10(exp(1))))
  
  # Append to list
  tables[[i]] <- model
}

# Combine all individual tables into one large dataframe
data <- do.call(rbind, tables)

# ---- Identify best models according to AIC and Maximum Likelihood ----

## Select best AIC (lowest AIC per subset and model)
best_models <- data %>%
  group_by(subset, modello) %>%
  slice_max(order_by = -AIC, n = 1)  # Note: slice_max with -AIC sorts ascending
  
colnames(best_models)[7] <- 'BEST_AIC'  # Rename column
best_models <- best_models[, c(3,4,7)]  # Keep only relevant columns

# Merge BEST_AIC info back into the full dataset
merged_df <- merge(data, best_models, by=c('modello', 'subset'))

## Select best Maximum Likelihood (highest MaxEstLhood)
best_models <- data %>%
  group_by(subset, modello) %>%
  slice_max(order_by = MaxEstLhood, n = 1)

colnames(best_models)[1] <- 'BEST_ML'  # Rename
best_models <- best_models[, c(1,3,4)]  # Keep relevant columns

# Merge BEST_ML info into the dataframe with BEST_AIC
merged_df <- merge(merged_df, best_models, by=c('modello', 'subset'))

# ---- Further subset: find overall best models ----

# Find the best model (minimum BEST_AIC) *per subset*
min_best_aic <- merged_df %>% 
  group_by(subset) %>% 
  slice_min(order_by = BEST_AIC, n = 1)

# Filter original dataframe to keep only those minimum BEST_AIC rows
min_best_aic_df <- merged_df %>% filter(BEST_AIC %in% min_best_aic$BEST_AIC)

# Find the absolute best model across all subsets
min_bestissimo_aic <- merged_df %>% 
  slice_min(order_by = BEST_AIC, n = 1)

min_bestissimo_aic_df <- merged_df %>% filter(BEST_AIC %in% min_bestissimo_aic$BEST_AIC)

# ---- Plot the AIC distributions ----

AIC <- ggplot(merged_df, aes(modello, AIC)) +
  geom_violin(width=0 ) +  # Draw violin plots without width
  facet_wrap(~subset, scales = 'free') +  # One facet per subset
  geom_point(aes(modello, BEST_AIC), color = 'red', size = 2) +  # Best AICs as red points
  geom_text(data = min_best_aic_df, aes(modello, BEST_AIC, label = "**"), color = 'black', size = 6, vjust = -4) +  # Mark subset-specific best models
  # geom_text(data = min_bestissimo_aic_df, aes(modello, BEST_AIC, label = "§"), color = 'black', size = 6, vjust = -3) + # Global best model - currently commented
  stat_summary(fun = mean, geom = "point", color = "blue", size = 2, shape=4) +  # Mean AIC as blue cross
  theme_bw() +
  scale_y_continuous(n.breaks = 20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        strip.background = element_rect(fill='white'), 
        strip.text = element_text(face = 'bold', size = 15)) +
  ylab('AIC') +
  xlab('') +
  ggtitle('Sao Miguel source')  # Title for plot

# Save the plot to a PDF
ggsave('AZM.pdf', AIC, 'pdf', width = 12, height = 7)

# ---- Final table with delta AIC and Akaike weights ----

# Get the best model per subset and model combination
grouped_data <- merged_df %>%
  group_by(subset, modello) %>%
  slice_max(order_by = -AIC, n = 1)

# Split into a list by subset
grouped_data <- grouped_data %>%
  group_by(subset) %>% 
  group_split()

# Process each subset separately to add deltaAIC, relative likelihoods and weights
final_list <- lapply(grouped_data, process_group)

# Combine into a final dataframe
final_df <- bind_rows(final_list)

# Save final results as a tab-separated file
write.table(final_df, "AZM-Final_table.tsv", quote = F, sep = "\t", row.names = F)
```

## STEP 6: Non-Parametric Bootstrapping (100 bootstraps) for parameter estimation on best selected model. 
Perform non-parametric bootstrapping to generate 100 observed SFS files by sampling SNPs from the SNP pool with replacement. 

The bootstrap procedure involves generating 100 SFS files and estimating parameters for the best model for each file.

For each bootstrap, adjust the `.tpl` file by adding the mutation rate in the last line.

Example of mutation rate adjustment in `.tpl` file:

```
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 2.1e-9  ## Set the mutation rate here
```

Create folders for each bootstrap with corresponding `.tpl`, `.est`, and `MSFS.obs` files, named according to the bootstrap iteration.

Folder organization and commands are as sthose of Step 5 but with a small difference within the fastsimcoal2 script:
```
 fsc27093 -t *.tpl -e *.est -m -n50000 -L30 -c10 -B10 -M --multiSFS
```

Once all bootstraps are run, we can calculate the 95% confidence intervals for the parameter estimates.

## STEP 7: Calculate 95% CI as in Marchi et al. (2023)  
```
# ---- Load required libraries ----
library(dplyr)      # For data manipulation
library(tidyr)      # For reshaping data (pivoting)

# ---- Set working directory ----
setwd('/home/user/parameter_estimation')

# ---- Load all files ending with 'bestlhoods' ----
samples_blh <- list.files(path = ".", 
                          pattern = "bestlhoods$", 
                          recursive = TRUE, 
                          full.names = TRUE)       

# Initialize empty list to store individual data frames
tables <- list()

# ---- Read and annotate each file ----
for (i in 1:length(samples_blh)) {
  model <- read.table(samples_blh[i], header = TRUE)  # Read the file into a dataframe
  
  # Extract 'bootstrap' ID from filename (e.g., bootstrap_1, bootstrap_2, etc.)
  boot_match <- regmatches(samples_blh[i], gregexec('bootstrap_[0-9]*', samples_blh[i]))
  model$bootstrap <- paste(unlist(boot_match)[1], collapse = " ")
  
  # Extract 'run' ID from filename (e.g., run1, run2, etc.)
  run_match <- regmatches(samples_blh[i], gregexec('run[0-9]*', samples_blh[i]))
  model$run <- paste(unlist(run_match)[1], collapse = " ")
  
  # Add the dataframe to the list
  tables[[i]] <- model
}

# ---- Merge all individual tables into a single dataframe ----
data <- do.call(rbind, tables)

# ---- Select the best run per bootstrap based on Maximum Estimated Likelihood ----
best_models <- data %>%
  group_by(bootstrap) %>%
  slice_max(order_by = MaxEstLhood, n = 1)  # Keep only the best run (highest MaxEstLhood) per bootstrap replicate

# ---- Select only specific parameter columns for further analysis ----
best_models_selected <- best_models[,c("N_POP0", "N_POP1", "N_POP2", "N_POP3", "N_POP4",
                                       "N_BOT1", "N_BOT2", "N_BOT3", "N_BOT4",
                                       "DURBOTUsa", "DURBOTAZm", "DURBOTAZj", "DURBOTIt", "DETTIME")]

# Define column names for easy reference
cols <- c("N_POP0", "N_POP1", "N_POP2", "N_POP3", "N_POP4",
          "N_BOT1", "N_BOT2", "N_BOT3", "N_BOT4",
          "DURBOTUsa", "DURBOTAZm", "DURBOTAZj", "DURBOTIt", "DETTIME")

# ---- Calculate bootstrap statistics: mean, 2.5% quantile, 97.5% quantile ----
ci_df <- best_models_selected %>%
  select(all_of(cols)) %>%                           # Select only target parameters
  pivot_longer(cols = everything(),                  # Reshape into long format (one row per parameter-value pair)
               names_to = "parameter", 
               values_to = "value") %>%
  group_by(parameter) %>%
  summarise(
    mean = mean(value),                               # Mean of bootstrap estimates
    q025 = quantile(value, 0.025),                    # 2.5% percentile (lower bound of 95% CI)
    q975 = quantile(value, 0.975),                    # 97.5% percentile (upper bound of 95% CI)
    .groups = "drop"
  )

# ---- Reshape again in long format (raw values) ----
# Useful for plotting distributions of bootstrap estimates
long_df <- best_models_selected %>%
  select(all_of(cols)) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")
```

### <u>CITATION</u>: 
Nina Marchi, Adamandia Kapopoulou, Laurent Excoffier. Demogenomic inference from spatially and temporally heterogeneous samples. *Molecular Ecology Resources*, 2023, 24 (1).
