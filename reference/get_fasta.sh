#!/bin/bash

#Get GRCh38 chr6
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr6.fa.gz
gunzip chr6.fa.gz

#Get NC_000006 which in chr6 from NCBI. Note this is build independent
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&amp;id=NC_000006&amp;rettype=fasta&amp;retmode=text" > NC_000006.fa

#Get GRCh38 transcript details
#wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
#gunzip refGene.txt.gz
