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



## TO DO
Remove warnings in python variant extraction script
Automate code for hg19 to hg38 liftover of gnomad r2_1 variants

