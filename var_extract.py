#!/usr/bin/python3

import os
import re
import pandas as pd
import numpy as np
from hgvs.easy import *

# Change to desired directory
os.chdir("/mnt/c/Users/hiphe/OneDrive - yes.my/Google Drive/Education/UCSI University/Final Year/Co-op/LekLab")
print("The current working directory is:", os.getcwd())

print('\
\n-----------------------------------\
\n\
\n         ClinVar Database          \
\n\
\n-----------------------------------')

# Load the *.vcf file
clin_var=pd.read_csv("./clinvar_20200905.vcf",
			sep="\t",
			skiprows=27,
			header=0)

# Filter out LAMA2 variants without conflicting interpretatation
lama2_clin_var = clin_var[clin_var.INFO.str.match(r'(.*LAMA2.*)')]
lama2_clin_var = lama2_clin_var[lama2_clin_var.INFO.str.match(r'(?!.*Conflicting_interpretations_of_pathogenicity.*)')]

# Extract important counts
sig = lama2_clin_var.INFO.str.extractall(r';(CLNSIG=.+?);') # variant significance
sig = sig.iloc[:,0].str.split(pat = "=", expand=True).reset_index(level=1, drop=True).iloc[:,1]
'''
#### What happened above? ####

iloc: subsets the dataframe by [row,col];
str.split: separates by 'pat = x' into separate columns;
reset_index: ensures new values can bind to 'lama2' DataFrame
iloc picks the important column

These steps are repeated in the next two code chunks below
-------------------------------
'''
type = lama2_clin_var.INFO.str.extractall(r';(CLNVC=.+?);') # variant type
type = type.iloc[:,0].str.split(pat = "=", expand=True).reset_index(level=1, drop=True).iloc[:,1]

mc = lama2_clin_var.INFO.str.extractall(r';(MC=.+?);') # molecular consequence
mc = mc.iloc[:,0].str.split(pat = "|", expand=True).reset_index(level=1, drop=True).iloc[:,1]

hgvs = lama2_clin_var.INFO.str.extractall(r';(CLNHGVS=.+?);') # genetic coordinates
hgvs = hgvs.iloc[:,0].str.split(pat = "=", expand=True).reset_index(level=1, drop=True).iloc[:,1]

# Binding important data points to original table
lama2_clin_var.loc[:,'sig'] = sig
lama2_clin_var.loc[:,'type'] = type
lama2_clin_var.loc[:,'mc'] = mc
lama2_clin_var.loc[:,'hgvs'] = hgvs

# Standardising data frame
lama2_clin_var = lama2_clin_var.rename(columns={'ID': 'clinvar'})
lama2_clin_var = lama2_clin_var.loc[:,['clinvar','hgvs','sig','type','mc']]

# Filtering pathogenic variants
pathogenic_clin_var = lama2_clin_var[lama2_clin_var.sig.str.match(r'.*[Pp]athogenic')]

# Printing variant counts

## by clinical significance
print('\nCounts by clinical significance\n')
tmp = lama2_clin_var.sig.value_counts(dropna=False)
print(tmp)
print('Total variants:', tmp.sum())

## by variant type
print('\nCounts by pathogenic variant type\n')
tmp = pathogenic_clin_var.type.value_counts(dropna=False)
print(tmp)
print('Total variants:', tmp.sum())

## by pathogenic molecular consequence
print('\nCounts by pathogenic molecular consequence\n')
tmp = pathogenic_clin_var.mc.value_counts(dropna=False)
print(tmp)
print('Total variants:', tmp.sum())

###---------------------------------------------------------------------------------------------------------

print('\
\n--------------------------------------\
\n\
\n             LOVD Database            \
\n\
\n--------------------------------------')

# Load the file
lovd = pd.read_csv(
	"./LOVD_full_download_LAMA2_2020-09-09_10.12.36.txt",
	header = 0,
	quotechar = '"',
	doublequote = True,
	sep = "\t",
	skiprows = 4110,
	nrows = 2055
	)

lovd.columns = lovd.columns.str.replace("['\{{','\}}']",'')
lovd.columns = lovd.columns.str.replace("VariantOnGenome/",'')
lovd = lovd.drop_duplicates(subset = 'DNA/hg38')
lovd['DNA/hg38'] = "NC_000006.12:" + lovd['DNA/hg38'].astype(str)

lovd = lovd.rename(
	columns = {
	'DBID':'lovd',
	'DNA/hg38':'hgvs',
	'ClinicalClassification':'sig',
	'type':'type'}
	)
lovd = lovd.loc[:,['lovd','hgvs','sig','type']]

# Filtering pathogenic variants
pathogenic_lovd = lovd[lovd.sig.str.match(r'.*[Pp]athogenic.*', na=False)]

# Print variant counts

## by clinical significance
print('\nCounts by clinical significance \n')
tmp = lovd.sig.value_counts(dropna=False)
print(tmp)
print('Total variants:', tmp.sum())

## by pathogenic variant type
print('\nCounts by pathogenic variant type\n')
tmp = pathogenic_lovd.type.value_counts(dropna=False)
print(tmp)
print('Total variants:', tmp.sum())

###-------------------------------------------------------------------------------------------

print('\
\n----------------------------------\
\n\
\n           EGL Database           \
\n\
\n----------------------------------')

# Load file
egl = pd.read_csv('./EmVClass.2020-Q3.csv', usecols = range(0,9), header = None)

# Rename columns accordingly
egl.columns = ['egl',
		'gene',
		'',
		'exon',
		'hgvs',
		'protein_change',
		'sig',
		'last reviewed',
		'alias listing']


# Filter LAMA2
lama2_egl = egl[egl.iloc[:,1].str.match(r'LAMA2')]
lama2_egl = lama2_egl.loc[:,['egl','hgvs','sig']]

lama2_egl.hgvs = lama2_egl.hgvs.apply(parse)

count = 0

for i in range(0,lama2_egl.shape[0]):
	try:
		lama2_egl.hgvs.iloc[i] = c_to_g(lama2_egl.hgvs.iloc[i])
	except:
		print('Exception:', lama2_egl.hgvs.iloc[i])
		count += 1
		print('Number of exceptions:', count)
		pass

lama2_egl.hgvs = lama2_egl.hgvs.apply(str)

# Pathogenic variants
pathogenic_egl = lama2_egl[lama2_egl.sig.str.match(r'[Pp]athogenic')]

# Print variant counts

## by clinical significance
print('\nCounts by clinical significance\n')
tmp = lama2_egl.sig.value_counts(dropna=False)
print(tmp)
print('Total variants:', tmp.sum())

#----------------------------------------------------------------------------------------------

print('\
\n-----------------------------------\
\n\
\n	Merge datasets together      \
\n\
\n-----------------------------------')

merged_list = lama2_clin_var.merge(lovd, on = 'hgvs', how = 'outer').merge(lama2_egl, on = 'hgvs', how = 'outer')
print('\nNumber of unique variants:', merged_list.hgvs.value_counts().sum())

merged_list = pathogenic_clin_var.merge(pathogenic_lovd, on = 'hgvs', how = 'outer').merge(pathogenic_egl, on = 'hgvs', how = 'outer')
print('\nNumber of unique pathogenic variants:', merged_list.hgvs.value_counts().sum())

overlap_all = pathogenic_clin_var.merge(pathogenic_lovd, on = 'hgvs', how = 'inner').merge(pathogenic_egl, on = 'hgvs', how = 'inner')
print('\nOverlapping variants in all datasets:', overlap_all.hgvs.value_counts().sum())

overlap_lovd_egl = pathogenic_lovd.merge(pathogenic_egl, on = 'hgvs', how = 'inner')
print('\nOverlapping variants in LOVD and EGL:', overlap_lovd_egl.hgvs.value_counts().sum())

overlap_lovd_clinvar = pathogenic_lovd.merge(pathogenic_clin_var, on = 'hgvs', how = 'inner')
print('\nOverlapping variants in ClinVar and LOVD:', overlap_lovd_clinvar.hgvs.value_counts().sum())

overlap_clinvar_egl = pathogenic_clin_var.merge(pathogenic_egl, on = 'hgvs', how = 'inner')
print('\nOverlapping variants in EGL and ClinVar:', overlap_clinvar_egl.hgvs.value_counts().sum())
