# LAMA2 related Congenital Muscular Dystrophy Prevelance


## Requirements
Python
```
pandas
pyhgvs
pyfaidx
BioPython
requests
```
R
```
dplyr
#install in R. install.packages("dplyr")
```

Data needed for pyhgvs
```
NCBI Genbank Chromosome 6 (NC_000006) 
GRCh38 Chromosome 6
Refseq transcript coordinates for GRCh38
```

## Installation using virtualenv
```
#Setup python environment
mkdir env
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt

#Get sequence files needed for pyhgvs
cd reference
./get_fasta.sh

```

## Mutation database sources
**ClinVar** - the 4th October 2020 data freeze was used  
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2020/clinvar_20201004.vcf.gz
  
**Leiden Open Variant Database (LOVD)** - Accessed in October 2020 and downloaded from the LAMA2 gene specific mutation database    
https://databases.lovd.nl/shared/download/all/gene/LAMA2
  
**Emory Genetics Laboratory (EGL)** - The 2020 Q4 data freeze was used  
https://www.egl-eurofins.com/emvclass/CSVArchives/EmVClass.2020-Q4.csv 

Downloaded files are included in this repository as input files to reproduce results. This is important as LOVD does not have data freezes.

```
data/clinvar_20200905.vcf.gz
data/EmVClass.2020-Q3.csv
data/LOVD_full_download_LAMA2_2020-09-09_10.12.36.txt
```

## Beta priors
The R script uses beta priors generated from the ExAC data set from our previous [publication](https://pubmed.ncbi.nlm.nih.gov/31105274)
and taken fron the github [repo](https://github.com/leklab/prevalence_estimation)
```
data/beta_parameter_prior_ExAC.txt
```

## Running the script
```
python3 var_extract.py
python3 get_gnomad_af.py > LAMA2_known_pathogenic_novel_lof_gnomad_af.tsv
./prevalence_estimate.r
```

## Output
```
Estimates including known pathogenic variants and novel gnomAD loss of function variants
ALL estimated prevalence (per million):  8.30285
ALL confidence interval with  95 % confidence:  6.27363 - 10.53518
NFE estimated prevalence (per million):  10.10345
NFE confidence interval with  95 % confidence:  6.741657 - 13.94393
AFR estimated prevalence (per million):  4.260902
AFR confidence interval with  95 % confidence:  2.052136 - 7.016912
AMR estimated prevalence (per million):  6.451412
AMR confidence interval with  95 % confidence:  3.447242 - 10.10158
EAS estimated prevalence (per million):  1.794676
EAS confidence interval with  95 % confidence:  0.6288023 - 3.360779


Conservative estimates including only known pathogenic variants
ALL estimated prevalence (per million):  2.542371
ALL confidence interval with  95 % confidence:  1.748451 - 3.44127
NFE estimated prevalence (per million):  3.204605
NFE confidence interval with  95 % confidence:  1.893433 - 4.755221
AFR estimated prevalence (per million):  0.9964744
AFR confidence interval with  95 % confidence:  0.3189813 - 1.924058
AMR estimated prevalence (per million):  3.026621
AMR confidence interval with  95 % confidence:  1.3644 - 5.134642
EAS estimated prevalence (per million):  0.38448
EAS confidence interval with  95 % confidence:  0.06072101 - 0.8871178
```

## TO DO
Remove warnings in python variant extraction script  
Automate code for hg19 to hg38 liftover of gnomad r2_1 variants

