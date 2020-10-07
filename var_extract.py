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
effect = lama2_clin_var.INFO.str.extractall(r';(CLNVC=.+?);') # variant type
effect = effect.iloc[:,0].str.split(pat = "=", expand=True).reset_index(level=1, drop=True).iloc[:,1]

mc = lama2_clin_var.INFO.str.extractall(r';(MC=.+?);') # molecular consequence
mc = mc.iloc[:,0].str.split(pat = "|", expand=True).reset_index(level=1, drop=True).iloc[:,1]

hgvs = lama2_clin_var.INFO.str.extractall(r';(CLNHGVS=.+?);') # genetic coordinates
hgvs = hgvs.iloc[:,0].str.split(pat = "=", expand=True).reset_index(level=1, drop=True).iloc[:,1]

# Binding important data points to original table
lama2_clin_var.loc[:,'sig'] = sig
lama2_clin_var.loc[:,'effect'] = effect
lama2_clin_var.loc[:,'mc'] = mc
lama2_clin_var.loc[:,'hgvs'] = hgvs

# Standardising data frame
lama2_clin_var = lama2_clin_var.rename(columns={'ID': 'clinvar'})
lama2_clin_var = lama2_clin_var.loc[:,['clinvar','hgvs','sig','effect','mc']]

# Filtering pathogenic variants
pathogenic_clin_var = lama2_clin_var[lama2_clin_var.sig.str.match(r'.*[Pp]athogenic')]

# Printing variant counts

## by clinical significance
print('\nCounts by clinical significance\n')
tmp = lama2_clin_var.sig.value_counts(dropna=False)
print(tmp)
print('Total variants:', tmp.sum())

## by variant type
print('\nCounts by pathogenic variant effect\n')
tmp = pathogenic_clin_var.effect.value_counts(dropna=False)
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
	'type':'effect'}
	)
lovd = lovd.loc[:,['lovd','hgvs','sig','effect']]

# Filtering pathogenic variants
pathogenic_lovd = lovd[lovd.sig.str.match(r'.*[Pp]athogenic.*', na=False)]

# Print variant counts

## by clinical significance
print('\nCounts by clinical significance \n')
tmp = lovd.sig.value_counts(dropna=False)
print(tmp)
print('Total variants:', tmp.sum())

## by pathogenic variant effect
print('\nCounts by pathogenic variant effect\n')
tmp = pathogenic_lovd.effect.value_counts(dropna=False)
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
EGL = pd.read_csv('./EmVClass.2020-Q3.csv', usecols = range(0,9), header = None)

# Rename columns accordingly
EGL.columns = ['egl',
		'gene',
		'',
		'exon',
		'hgvs',
		'protein_change',
		'sig',
		'last reviewed',
		'alias listing']


# Filter LAMA2
lama2_EGL = EGL[EGL.iloc[:,1].str.match(r'LAMA2')]
lama2_EGL = lama2_EGL.loc[:,['egl','hgvs','sig']]

lama2_EGL.hgvs = lama2_EGL.hgvs.apply(parse)

try:
	for i in range(0,lama2_EGL.shape[0]):
		lama2_EGL.hgvs.iloc[i] = str(c_to_g(lama2_EGL.hgvs.iloc[i]))
except:
	pass
	#add a count line and print number of exceptions

# Pathogenic variants
pathogenic_EGL = lama2_EGL[lama2_EGL.sig.str.match(r'[Pp]athogenic')]

# Print variant counts

## by clinical significance
print('\nCounts by clinical significance\n')
tmp = lama2_EGL.sig.value_counts(dropna=False)
print(tmp)
print('Total variants:', tmp.sum())

#----------------------------------------------------------------------------------------------

print('\
\n-----------------------------------\
\n\
\n	Merge datasets together      \
\n\
\n-----------------------------------')
