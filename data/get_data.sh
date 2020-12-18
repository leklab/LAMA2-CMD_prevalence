#!/bin/bash

wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20200905.vcf.gz
curl https://databases.lovd.nl/shared/download/all/gene/LAMA2 > LOVD_full_download_LAMA2.txt
wget https://www.egl-eurofins.com/emvclass/CSVArchives/EmVClass.2020-Q4.csv
wget https://raw.githubusercontent.com/leklab/prevalence_estimation/master/data/beta_parameter_prior_ExAC.txt
