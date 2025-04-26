# VARIANT CALLING

The following document describe the procedure from raw WGS reads to the two final snp datasets used in population analyses.
The procedure follows the tutorial by Mark Ravinet & Joana Meier associated to the Physalia course on 'Speciation & Population Genomics'. The latter can be found at https://speciationgenomics.github.io/



## TRIMMING

The following commands will perform quality trimming on the raw Illumina reads (.fastq.gz). R1 and R2 reads for one individual (i.e. the raw data) are here conventionally named XXX_1.fastq.gz and XXX_2.fastq.gz. The commands are to be repeated on all individuals.

### Required software
fastp

### Input files
`XXX_1.fastq.gz` Raw data, R1 reads
`XXX_2.fastq.gz` Raw data, R2 reads

### Commands
```
# trimming
fastp \
--thread 10 \
-i XXX_1.fastq.gz -I XXX_2.fastq.gz \
-o XXX_1t.fastq.gz -O XXX_2t.fastq.gz \
-l 50 \
-h XXX_stats \
--dont_eval_duplication \
--correction \
--cut_right \
--cut_right_window_size 4 \
--cut_right_mean_quality 24 \
--cut_front \
--cut_front_window_size 1 \
--cut_front_mean_quality 20 \
--cut_tail \
--cut_tail_window_size 1 \
--cut_tail_mean_quality 20
```

### Relevant output files
`XXX_1t.fastq.gz` Trimmed data, R1 reads
`XXX_2t.fastq.gz` Trimmed data, R2 reads
`XXX_stats` Read quality and trimming statistics






## REMAPPING AND READ FILTERING

Reads are here remapped over the genome. The commands are to be repeated on all individuals.

### Required software
bbmap
samtools
picard

### Input files
`Masked.fasta` The genome sequence.
`Masked.fasta.fai` Its index.
`XXX_1.fastq.gz` Trimmed data, R1 reads
`XXX_2.fastq.gz` Trimmed data, R2 reads

### Commands
```
# remapping
bbmap.sh \
ref=Masked.fasta \
in=XXX_1t.fastq.gz \
in2=$XXX_2t.fastq.gz \
outm=XXX.mapped.bam \
t=20 \
maxindel=200 \
pairlen=500 \
1>XXX_out.txt 2>XXX_err.txt

# sort, index
samtools sort -o XXX.mapped.sort.bam XXX.mapped.bam 
samtools index -b XXX.mapped.sort.bam

# remove duplicate reads
java -jar picard.jar \
MarkDuplicates \
I=XXX.mapped.sort.bam \
O=XXX.mapped.sort.rmd.bam \
M=XXX.dupstats \
VALIDATION_STRINGENCY=SILENT \
REMOVE_DUPLICATES=true

# remove imperfectly matching pairs
samtools view -@ 10 -m 16G \
-q 20 \
-f 0x0002 \
-F 0x0004 \
-F 0x0008 \
-b XXX.mapped.sort.rmd.bam \
> XXX.mapped.sort.rmd.q20Ff.bam

# index
samtools index -b XXX.mapped.sort.rmd.q20Ff.bam

# calculate coverage statistics
pileup.sh \
in=XXX.mapped.sort.rmd.q20Ff.bam \
2>XXX_covstats.txt
```

### Relevant output files
`XXX.mapped.sort.rmd.q20Ff.bam` Reads of one individual mapped over genome, deprived of duplicate reads and imperfectly matching read pairs.






## GLOBAL CALLING

Individual remapping files are globally compared. All variable sites are initially identified and snp as well as indel variants are called using the multiallelic caller. Snps within the range of 3bp from an indel are removed. Indels are removed, and only snp retained for further analysis.

### Required software
bcftools

### Input files
`Masked.fasta` The genome sequence.
`Masked.fasta.fai` Its index.
`bamlist83` List of .bam files from previous step (e.g. XXX.mapped.sort.rmd.q20Ff.bam), one per individual

### Commands
```
# calling
bcftools mpileup \
-Ou \
-d 200 \
-a AD,DP,SP \
-f Masked.fasta \
--bam-list bamlist83 \
| bcftools call \
-f GQ,GP \
-Ou \
-m -v \
--ploidy 2 \
| bcftools filter \
--SnpGap 3:indel,other \
-Ou \
| bcftools view \
--threads 4 \
-v snps \
-Oz \
-o allsnps.vcf.gz

# trick!
# option -R in mpileup can be used to inlcude only a subset of contigs for trial runs/troubleshooting

# sort, index
bcftools sort allsnps.vcf.gz -Oz -o allsnps.sort.vcf.gz
bcftools index allsnps.sort.vcf.gz

# explore
bcftools stats allsnps.sort.vcf.gz > allsnps.sort.vcf.gz_stats
```

### Relevant output files
`allsnps.sort.vcf.gz` Called snps, to be filtered.





## FILTERING FOR CONTIG/REGIONS

Initial filtering is conducted on the snp set to remove a) sites in regions identified as repeated elements and b) sites in contigs characterized by an extreme coverage.

### Required software
vcftools
bcftools

### Input files
`allsnps.sort.vcf.gz` Called snps (to be filtered).
`RepeatElements.bed` List of regions identified as repeated elements (to be excluded).
`contigbuoni_inner` Contigs characterized by a standard coverage (to be reatained).

### Commands
```
# exclude snps in repeated elements
vcftools --gzvcf allsnps.sort.vcf.gz \
--exclude-bed RepeatElements.bed \
--recode \
--stdout | gzip -c > allsnps.sort.nomask.vcf.gz

# index
bcftools index allsnps.sort.nomask.vcf.gz

# exclude snps in contigs with extreme coverage
bcftools view \
allsnps.sort.nomask.vcf.gz \
--threads 12 \
-R contigbuoni_inner \
-Oz \
-o allsnps.sort.nomask.nobadcov.vcf.gz

# index
bcftools index allsnps.sort.nomask.nobadcov.vcf.gz

# explore
bcftools stats allsnps.sort.nomask.nobadcov.vcf.gz > allsnps.sort.nomask.nobadcov.vcf.gz_stats
```

### Relevant output files
`allsnps.sort.nomask.nobadcov.vcf.gz` Called snps, filtered for contigs/regions.




## QUALITY FILTERING OF SITES

The distribution of multiple statistics are calculated over sites and over individuals to identify sensible filtering parameters/thresholds. These filters are then applied to the snp list from the previous step to produce a final, filtered, dataset. Linkage pruning is applied to the latter to produce a twin dataset with sites in linkage equilibrium, to be used in some analyses.

### Required software
vcftools
bcftools
plink

### Input files
`allsnps.sort.nomask.nobadcov.vcf.gz`

### Commands
```
# calculate the distribution of multiple snp statistics
# manually revise to identify sensible filtering thresholds
# allele frequency
vcftools --gzvcf allsnps.sort.nomask.nobadcov.vcf.gz --freq2 --max-alleles 2 --stdout > snps.allelefreq
# mean depth per individual
vcftools --gzvcf allsnps.sort.nomask.nobadcov.vcf.gz --depth --stdout > snps.MeanDepthXind
# mean depth per site
vcftools --gzvcf allsnps.sort.nomask.nobadcov.vcf.gz --site-mean-depth --stdout > snps.sitemeandepth
# site quality
vcftools --gzvcf allsnps.sort.nomask.nobadcov.vcf.gz --site-quality --stdout > snps.sitequality
# missing sites per individual
vcftools --gzvcf allsnps.sort.nomask.nobadcov.vcf.gz --missing-indv --stdout > snps.missingSitesxInd
# missing individuals per site
vcftools --gzvcf allsnps.sort.nomask.nobadcov.vcf.gz --missing-site --stdout > snps.missingIndxSite
# heterozygosity
vcftools --gzvcf allsnps.sort.nomask.nobadcov.vcf.gz --het --stdout > snps.heterozygosity

# apply filters
vcftools --gzvcf allsnps.sort.nomask.nobadcov.vcf.gz \
--remove-indels \
--thin 5 \
--min-alleles 2 \
--max-alleles 2 \
--maf 0.015 \
--max-missing 0.95 \
--minQ 50 \
--min-meanDP 15 \
--max-meanDP 52 \
--minDP 15 \
--maxDP 52 \
--recode \
--stdout | gzip -c > snps_3p.vcf.gz

# index
bcftools index snps_3p.vcf.gz

# explore 
bcftools stats snps_3p.vcf.gz > snps_3p.vcf.gz_stats

# identify unlinked sites
plink \
--vcf snps_3p.vcf.gz \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out snps_3p

# remove unlinked sites
# at bonus, it produces input files for pca and admixture
plink \
--vcf snps_3p.vcf.gz \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--extract snps_3p.prune.in \
--make-bed \
--pca \
--out snps_3p_unl

# recode snps as uncompressed vcf
plink \
--bfile snps_3p_unl \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--recode vcf \
--out snps_3p_unl

# compress vcf
bcftools view -I snps_3p_unl.vcf -Oz -o snps_3p_unl.vcf.gz

# index
bcftools index snps_3p_unl.vcf.gz

# explore 
bcftools stats snps_3p_unl.vcf.gz > snps_3p_unl.vcf.gz_stats
```

### Relevant output files
`snps_3p.vcf.gz` Called snps, final filtered version, all sites.
`snps_3p_unl.vcf.gz` Called snps, final filtered version, unlinked sites only.


