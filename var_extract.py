#!/usr/bin/python3
'''
Necessary packages/modules:
- pandas
- numpy
- hgvs (https://github.com/biocommons/hgvs)

Dependencies:
- ClinVar data freeze (weekly, *.vcf)
- LAMA2 LOVD data freeze (quarterly, *.txt)
- EGL data freeze (quarterly, *.csv)
- For pyhgvs functionality: refGene file from UCSC (ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz) [edit the LAMA2 transcript name to NM_000426.3]
- For pyhgvs functionality: GRCh38 from both NCBI and UCSC (to conserve space, download just chr6)
'''
# Load necessary modules
import os
import re
import pandas as pd
import numpy as np
import pyhgvs
from pyfaidx import Fasta
import pyhgvs.utils as hgvs_utils

# Change to desired directory
os.chdir("/mnt/c/Users/hiphe/OneDrive - yes.my/Google Drive/Education/UCSI University/Final Year/Co-op/LekLab")
print("The current working directory is:", os.getcwd())

# Important dictionaries
dsig = {'VOUS':'VUS',
	'Uncertain_significance':'VUS'}

dtype = {'single_nucleotide_variant':'SNV',
	'subst':'SNV',
	'Duplication':'dup',
	'Indel':'delins',
	'Insertion':'ins',
	'Deletion':'del',
	'Duplication':'dup'}

# Fetching reference files for HGVS to VCF genomic coordinate formatting
genome = Fasta('/mnt/c/Users/hiphe/Downloads/GRCh38_latest_genomic.fna')
genome2 = Fasta('/mnt/c/Users/hiphe/Downloads/chr6.fa')
with open('./refGene.txt') as infile:
	transcripts = hgvs_utils.read_transcripts(infile)

def get_transcript(name):
	return transcripts.get(name)

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

# Filter for LAMA2
lama2_clin_var = clin_var[clin_var.INFO.str.match(r'(.*LAMA2.*)')]
lama2_clin_var['ID'] = 'clinvar_' + lama2_clin_var['ID'].astype(str)

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
lama2_clin_var = lama2_clin_var.loc[:,['ID','CHROM','POS','REF','ALT','hgvs','sig','type','mc']]
lama2_clin_var = lama2_clin_var.replace({'type':dtype})
lama2_clin_var = lama2_clin_var.replace('Conflicting_interpretations_of_pathogenicity','VUS')
lama2_clin_var = lama2_clin_var.replace({'sig':dsig})

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
lovd_ori = pd.read_csv(
	"./LOVD_full_download_LAMA2_2020-09-09_10.12.36.txt",
	header = 0,
	quotechar = '"',
	doublequote = True,
	sep = "\t",
	skiprows = 4110,
	nrows = 2055
	)

# Cleaning up the names of columns
lovd_ori.columns = lovd_ori.columns.str.replace("['\{{','\}}']",'')
lovd_ori.columns = lovd_ori.columns.str.replace("VariantOnGenome/",'')
lovd_ori = lovd_ori.drop_duplicates(subset = 'DNA/hg38')
lovd_ori['DNA/hg38'] = "NC_000006.12:" + lovd_ori['DNA/hg38'].astype(str)

lovd = lovd_ori.rename(
		columns = {
		'DBID':'ID',
		'DNA/hg38':'hgvs',
		'ClinicalClassification':'sig'}
		)

# Selecting important columns
lovd = lovd.loc[:,['ID','hgvs','sig','type']]
lovd = lovd.replace({'type':dtype})

# Converting HGVS to VCF genomic coordinates
count = 0

for i in lovd.index:
	try:
		CHROM, POS, REF, ALT = pyhgvs.parse_hgvs_name(lovd.loc[i,'hgvs'], genome)
	except:
		print('\nException:', lovd.loc[i,'hgvs'])
		count += 1
		print('Number of exceptions:', count)
		pass

	lovd.loc[i,'CHROM'] = 6
	lovd.loc[i,"POS"] = POS
	lovd.loc[i,"REF"] = REF
	lovd.loc[i,"ALT"] = ALT

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
egl.columns = ['ID',
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
lama2_egl = lama2_egl.loc[:,['ID','hgvs','sig']]
lama2_egl = lama2_egl.replace({'sig':dsig})
lama2_egl['ID'] = 'egl_' + lama2_egl['ID'].astype(str)

# Converting HGVS to VCF genomic coordinates
count = 0

for i in lama2_egl.index:
	try:
		CHROM, POS, REF, ALT = pyhgvs.parse_hgvs_name(lama2_egl.loc[i,'hgvs'],
								genome2,
								get_transcript=get_transcript)
	except:
		print('\nException:', lama2_egl.loc[i,'hgvs'])
		count += 1
		print('Number of exceptions:', count)
		pass

	lama2_egl.loc[i,"CHROM"] = 6
	lama2_egl.loc[i,"POS"] = POS
	lama2_egl.loc[i,"REF"] = REF
	lama2_egl.loc[i,"ALT"] = ALT

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

# Full list of LAMA2 variants
all_vars = lama2_clin_var.merge(lovd, on = ["CHROM","POS","REF","ALT"], how = 'outer').merge(lama2_egl, on = ["CHROM","POS","REF","ALT"], how = 'outer')

## Resolving conflicts among significances
for i in all_vars.index :

	sig_clinvar = all_vars.loc[i,'sig_x']
	sig_lovd = all_vars.loc[i,'sig_y']
	sig_egl = all_vars.loc[i,'sig']
	sigs = [keep for keep in [sig_clinvar, sig_lovd, sig_egl] if str(keep) != 'nan'] #keep only the non-NaN values

	if ('VUS' in sigs) :
		all_vars.loc[i,'sig'] = 'VUS'

	elif all(re.match(r'.*([Bb]enign).*', sig) for sig in sigs) :
		all_vars.loc[i,'sig'] = 'Benign'

	elif all(re.match(r'.*([Pp]athogenic).*', sig) for sig in sigs) :
		all_vars.loc[i,'sig'] = 'Pathogenic'

	elif np.unique(sig).size != len(sig) : #shows that each element sig is unique
		all_vars.loc[i,'sig'] = 'VUS'
	else:
		print("There's an exception here:\n",
			"Index:", i, '\n',
			'Issue resolving: sigs \n')

## Display value counts for variant significance
print(all_vars.sig.value_counts(), '\nTotal variants:', all_vars.shape[0])

## Standardizing columns to VCF format
all_vars['ID'] = all_vars.loc[:,['ID_x','ID_y','ID']].apply(lambda x:
								','.join(x.dropna().astype(str)), axis=1)
all_vars['QUAL'] = '.'
all_vars['FILTER'] = '.'
all_vars['INFO'] = all_vars.loc[:,
				['sig','type_x','type_y']
				 ].apply(
				  	lambda x: ','.join(x.dropna().astype(str)), axis=1
				  	)

all_vars = all_vars[['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']].sort_values(by = 'POS')

# Pathogenic only
pathogenic_all = all_vars[all_vars.INFO.str.match(r'Pathogenic')]
print('\nNumber of unique pathogenic variants:', pathogenic_all.shape[0])

# Overlapping between all
overlap_all = pathogenic_clin_var.merge(pathogenic_lovd, on = ["CHROM","POS","REF","ALT"], how = 'inner').merge(pathogenic_egl, on = ["CHROM","POS","REF","ALT"], how = 'inner')
print('\nOverlapping variants in all datasets:', overlap_all.shape[0])

# Overlapping between LOVD and EGL
overlap_lovd_egl = pathogenic_lovd.merge(pathogenic_egl, on = ["CHROM","POS","REF","ALT"], how = 'inner')
print('\nOverlapping variants in LOVD and EGL:', overlap_lovd_egl.shape[0])

# Overlapping between LOVD and ClinVar
overlap_lovd_clinvar = pathogenic_lovd.merge(pathogenic_clin_var, on = ["CHROM","POS","REF","ALT"], how = 'inner')
print('\nOverlapping variants in ClinVar and LOVD:', overlap_lovd_clinvar.shape[0])

# Overlapping between ClinVar and EGL
overlap_clinvar_egl = pathogenic_clin_var.merge(pathogenic_egl, on = ["CHROM","POS","REF","ALT"], how = 'inner')
print('\nOverlapping variants in EGL and ClinVar:', overlap_clinvar_egl.shape[0])
