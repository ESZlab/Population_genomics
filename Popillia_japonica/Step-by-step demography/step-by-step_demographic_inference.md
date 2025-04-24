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

Run 50 independent runs for each model.   ***********************QUI DA SISTEMARE CON I LOOP A 50 volte 

```
fsc27093 -t model*.tpl -e model*.est -m -n500000 -L40 -c12 -B12 --removeZeroSFS -M --multiSFS
```

### SUGGESTION: 
Quickly run the code with `-L 1` and `-n10000` to generate a `.par` file for model verification.

Then, convert the `.par` file to a `.pdf` using the script available at: https://cmpg.unibe.ch/software/fastsimcoal26/parFileInterpreter.html.

## STEP 5: Model Selection and Parameter Estimation
After running Fastsimcoal2, we analyze the results and select the best model according to model comparison with the Akaike Information Criterion (AIC). **************QUI ANDREBBE AGGIUNTO QUELLO CHE SI é FATTO e lo script di R

## STEP 6: Non-Parametric Bootstrapping (100 bootstraps) for parameter estimation on best selected model. 
Perform non-parametric bootstrapping to generate 100 observed SFS files by sampling SNPs from the SNP pool with replacement.  *****************NON HO LO SCRIPT DIETRO, DOVREBBE ESSERE NEL MIO PC NEL LAB, SPERO DI AVERLO SCRITTO O CHE SIA NELL'HISTORY, CERCANDO 'SHUF', 'BOOTSTRAP'. Anche se forse non serve metterlo

The bootstrap procedure involves generating 100 SFS files and estimating parameters for the best model for each file.

For each bootstrap, adjust the `.tpl` file by adding the mutation rate in the last line.

Example of mutation rate adjustment in `.tpl` file:

```
//per block: Datatype, numm loci, rec rate and mut rate + optional parameters
FREQ 1 0 2.1e-9  ## Set the mutation rate here
```

Create folders for each bootstrap with corresponding `.tpl`, `.est`, and `MSFS.obs` files, named according to the bootstrap iteration.
# Command line    *****************************qui ci va il comando di fastsimcoal2 con il loop per 100 bootstrap

Once all bootstraps are run, we can calculate the 95% confidence intervals for the parameter estimates.

## STEP 7: Calculate 95% CI as in Marchi et al. (2023)    ********************** qui ci va lo script?

### <u>CITATION</u>: 
Nina Marchi, Adamandia Kapopoulou, Laurent Excoffier. Demogenomic inference from spatially and temporally heterogeneous samples. *Molecular Ecology Resources*, 2023, 24 (1).
