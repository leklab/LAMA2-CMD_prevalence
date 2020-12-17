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
import subprocess

# Change to desired directory
#os.chdir("/mnt/c/Users/hiphe/OneDrive - yes.my/Google Drive/Education/UCSI University/Final Year/Co-op/LekLab")
#print("The current working directory is:", os.getcwd())

'''
-----------------------------------

-   	 Variant Aggregation      -

-----------------------------------
'''

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
'''
genome = Fasta('/mnt/c/Users/hiphe/Downloads/GRCh38_latest_genomic.fna')
genome2 = Fasta('/mnt/c/Users/hiphe/Downloads/chr6.fa')
'''

with open('reference/refGene.txt') as infile:
	transcripts = hgvs_utils.read_transcripts(infile)

def get_transcript(name):
	return transcripts.get(name)


def clinvar_variants(filename) :
	print('\
	\n-----------------------------------\
	\n\
	\n         ClinVar Database          \
	\n\
	\n-----------------------------------')

	# Load the *.vcf.gz file
	clinvar=pd.read_csv(filename,
				compression='gzip',
				sep="\t",
				skiprows=27,
				header=0)

	# Filter for LAMA2
	lama2_clinvar = clinvar[clinvar.INFO.str.match(r'(.*LAMA2.*)')]
	lama2_clinvar['ID'] = 'clinvar_' + lama2_clinvar['ID'].astype(str)

	# Extract important counts
	sig = lama2_clinvar.INFO.str.extractall(r';(CLNSIG=.+?);') # variant significance
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
	type = lama2_clinvar.INFO.str.extractall(r';(CLNVC=.+?);') # variant type
	type = type.iloc[:,0].str.split(pat = "=", expand=True).reset_index(level=1, drop=True).iloc[:,1]

	mc = lama2_clinvar.INFO.str.extractall(r';(MC=.+?);') # molecular consequence
	mc = mc.iloc[:,0].str.split(pat = "|", expand=True).reset_index(level=1, drop=True).iloc[:,1]

	hgvs = lama2_clinvar.INFO.str.extractall(r';(CLNHGVS=.+?);') # genetic coordinates
	hgvs = hgvs.iloc[:,0].str.split(pat = "=", expand=True).reset_index(level=1, drop=True).iloc[:,1]

	# Binding important data points to original table
	lama2_clinvar.loc[:,'sig'] = sig
	lama2_clinvar.loc[:,'type'] = type
	lama2_clinvar.loc[:,'mc'] = mc
	lama2_clinvar.loc[:,'hgvs'] = hgvs

	#Remove # in front of CHROM as column label
	lama2_clinvar.rename(columns={'#CHROM': 'CHROM'},inplace=True)

	# Standardising data frame
	#lama2_clinvar = lama2_clinvar.loc[:,['ID','CHROM','POS','REF','ALT','hgvs','sig','type','mc']]
	lama2_clinvar = lama2_clinvar.reindex(columns=['ID','CHROM','POS','REF','ALT','hgvs','sig','type','mc'])
	lama2_clinvar = lama2_clinvar.replace({'type':dtype})
	lama2_clinvar = lama2_clinvar.replace('Conflicting_interpretations_of_pathogenicity','VUS') # all conflict resolved to 'uncertain significance'
	lama2_clinvar = lama2_clinvar.replace({'sig':dsig})

	# Filtering pathogenic variants
	#pathogenic_clinvar = lama2_clinvar[lama2_clinvar.sig.str.match(r'.*[Pp]athogenic')]
	#return [lama2_clinvar, pathogenic_clinvar]
	#print(lama2_clinvar)

	return lama2_clinvar
###---------------------------------------------------------------------------------------------------------


def lovd_variants(filename) :

	print('\
	\n--------------------------------------\
	\n\
	\n             LOVD Database            \
	\n\
	\n--------------------------------------')

	# Load the file
	lovd_ori = pd.read_csv(
		filename,
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


	genome = Fasta('reference/NC_000006.fa')

	# Converting HGVS to VCF genomic coordinates
	count = 0

		
	for i in lovd.index:
		try:
			#print('Trying %d\n' %(lovd.loc[i,'hgvs']))
			CHROM, POS, REF, ALT = pyhgvs.parse_hgvs_name(lovd.loc[i,'hgvs'], genome)
			lovd.loc[i,'CHROM'] = 6
			lovd.loc[i,"POS"] = POS
			lovd.loc[i,"REF"] = REF
			lovd.loc[i,"ALT"] = ALT
		except:
			print('\nException:', lovd.loc[i,'hgvs'])
			count += 1
			print('Number of exceptions:', count)
			pass
	
	# Filtering pathogenic variants
	#pathogenic_lovd = lovd[lovd.sig.str.match(r'.*[Pp]athogenic.*', na=False)]

	#return [lovd, pathogenic_lovd]
	lovd.dropna(inplace=True)
	lovd = lovd.astype({'CHROM': 'int32', 'POS': 'int32'})

	#print(lovd)

	return lovd
###-------------------------------------------------------------------------------------------



def egl_variants(filename) :
	print('\
	\n----------------------------------\
	\n\
	\n           EGL Database           \
	\n\
	\n----------------------------------')

	# Load file
	egl = pd.read_csv(filename, usecols = range(0,9), header = None)

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

	genome = Fasta('reference/chr6.fa')

	# Converting HGVS to VCF genomic coordinates
	count = 0
	transcript = get_transcript('NM_000426')

	
	for i in lama2_egl.index:
		try:
			#CHROM, POS, REF, ALT = pyhgvs.parse_hgvs_name(lama2_egl.loc[i,'hgvs'], genome, get_transcript=get_transcript)
			CHROM, POS, REF, ALT = pyhgvs.parse_hgvs_name(lama2_egl.loc[i,'hgvs'], genome, transcript)
			lama2_egl.loc[i,"CHROM"] = 6
			lama2_egl.loc[i,"POS"] = POS
			lama2_egl.loc[i,"REF"] = REF
			lama2_egl.loc[i,"ALT"] = ALT

		except:
			print('\nException:', lama2_egl.loc[i,'hgvs'])
			count += 1
			print('Number of exceptions:', count)
			pass
	
	lama2_egl = lama2_egl.astype({'CHROM': 'int32', 'POS': 'int32'})

	# Pathogenic variants
	#pathogenic_egl = lama2_egl[lama2_egl.sig.str.match(r'[Pp]athogenic')]

	#return [lama2_egl, pathogenic_egl]
	return lama2_egl
#----------------------------------------------------------------------------------------------


'''
def counts(database, query) :

	# Determine database file location
	if database == 'clinvar':
		file = input('Provide ClinVar database location: ')
		dataset = clinvar_variants(file)

	elif database == 'lovd':
		file = input('Provide LOVD database location: ')
		dataset = lovd_variants(file)

	elif database == 'egl':
		file = input('Provide EGL database location: ')
		dataset = egl_variants(file)

	else:
		print('Compatible databases: clinvar, lovd, egl')

	# Determine which variant count to give
	if query == 'significance':
		significance = dataset[0]
		print('\nCounts by clinical significance\n')
		tmp = significance.sig.value_counts(dropna=False)
		print(tmp)
		print('Total variants:', tmp.sum())

	elif query == 'type':
		vtype = dataset[1]
		print('\nCounts by pathogenic variant type\n')
		tmp = vtype.type.value_counts(dropna=False)
		print(tmp)
		print('Total variants:', tmp.sum())

	else:
		print('\nCompatible queries: significance, type')

#-----------------------------------------------------------------------------------------------
'''


def merge_datasets(clinvar, lovd, egl, output):
	print('\
	\n-----------------------------------\
	\n\
	\n	Merge datasets together      \
	\n\
	\n-----------------------------------')

	# Full list of LAMA2 variants
	all_vars = clinvar.merge(lovd, on = ["CHROM","POS","REF","ALT"], how = 'outer').merge(egl, on = ["CHROM","POS","REF","ALT"], how = 'outer')

	## Resolving conflicts among significances
	for i in all_vars.index :

		sig_clinvar = all_vars.loc[i,'sig_x']
		sig_lovd = all_vars.loc[i,'sig_y']
		sig_egl = all_vars.loc[i,'sig']
		sigs = [keep for keep in [sig_clinvar, sig_lovd, sig_egl] if str(keep) != 'nan'] #keep only the non-NaN values

		if ('VUS' in sigs) :
			all_vars.loc[i,'sig'] = 'VUS'

		elif all(re.match(r'.*([Bb]enign).*', val) for val in sigs) :
			all_vars.loc[i,'sig'] = 'Benign'

		elif all(re.match(r'.*([Pp]athogenic).*', val) for val in sigs) :
			all_vars.loc[i,'sig'] = 'Pathogenic'

		elif np.unique(sigs).size != len(sigs) : #shows that each element sig is unique
			all_vars.loc[i,'sig'] = 'VUS'

		else :
			print("There's an exception here:\n",
				"Index:", i, '\n',
				'Issue resolving:', sigs ,'\n')

	## Display value counts for variant significance
	print(all_vars.sig.value_counts(), '\nTotal variants:', all_vars.shape[0])

	## Standardizing columns to VCF format
	all_vars['ID'] = all_vars[['ID_x','ID_y','ID']].apply(lambda x: '|'.join(x.dropna().astype(str)), axis=1)
	all_vars['QUAL'] = '.'
	all_vars['FILTER'] = '.'
	all_vars['INFO'] = all_vars[['sig','type_x','type_y']].apply(lambda x: '|'.join(x.dropna().astype(str)), axis=1)

	all_vars = all_vars[['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']].sort_values(by = 'POS')

	all_vars.to_csv('LAMA2_unique_variants.tsv', sep = '\t', header = False, index = False)

	# Pathogenic only
	pathogenic_all = all_vars[all_vars.INFO.str.match(r'Pathogenic')]
	print('\nNumber of unique pathogenic variants:', pathogenic_all.shape[0])

	# Overlapping variant counts

	# In all datasets
	overlap_all = pathogenic_all[(pathogenic_all.ID.str.match(r'.*(clinvar).*(LAMA2).*(egl).*'))]
	print('Overlapping variants in all datasets:', overlap_all.shape[0])

	# In LOVD and EGL
	overlap_lovd_egl = pathogenic_all[pathogenic_all.ID.str.match('LAMA2_\d+\|.*')]
	print('Overlapping variants in LOVD and EGL:', overlap_lovd_egl.shape[0])

	# In ClinVar and LOVD
	overlap_clinvar_lovd = pathogenic_all[pathogenic_all.ID.str.match('clinvar_\d+\|LAMA2_\d+$')]
	print('Overlapping variants in ClinVar and LOVD:', overlap_clinvar_lovd.shape[0])

	# In ClinVar and EGL
	overlap_clinvar_egl = pathogenic_all[pathogenic_all.ID.str.match(r'clinvar_\d+\|(?!LAMA2_\d+)')]
	print('Overlapping variants in EGL and ClinVar:', overlap_clinvar_egl.shape[0])

	# In ClinVar only
	clinvar_only = pathogenic_all[pathogenic_all.ID.str.match('clinvar_\d+$')]
	print('ClinVar-exclusive variants:', clinvar_only.shape[0])

	# In LOVD only
	lovd_only = pathogenic_all[pathogenic_all.ID.str.match('LAMA2_\d+$')]
	print('LOVD-exclusive variants:', lovd_only.shape[0])

	# In EGL only
	egl_only = pathogenic_all[pathogenic_all.ID.str.match('egl_\d+')]
	print('EGL-exclusive variants:', egl_only.shape[0])

	# Checking that all variants counted
	sum_total = overlap_all.shape[0]\
			+ overlap_clinvar_lovd.shape[0]\
			+ overlap_lovd_egl.shape[0]\
			+ overlap_clinvar_egl.shape[0]\
			+ clinvar_only.shape[0]\
			+ lovd_only.shape[0]\
			+ egl_only.shape[0]\

	print('\nThe True Test: variant counts all good?')
	if sum_total == pathogenic_all.shape[0]:
		print('\nAll variants accounted for')
	else:
		print('\nBetter check your code, bro')

	pathogenic_all.to_csv(output, sep = '\t', header = False, index = False)
#-------------------------------------------------------------------


'''
print('\
\n-----------------------------------\
\n\
\n	 Prevalence Estimation       \
\n\
\n-----------------------------------')

'''

# Estimate without novel gnomAD alleles
#subprocess.check_call(['Rscript', './prevalence_estimate.r'], shell = False)

if __name__ == '__main__' :

	clinvar = clinvar_variants('data/clinvar_20200905.vcf.gz')
	#clinvar.to_csv('test.tsv',index=False,sep='\t')

	lovd = lovd_variants(filename = 'data/LOVD_full_download_LAMA2_2020-09-09_10.12.36.txt')
	#lovd.to_csv('lovd2.tsv',index=False,sep='\t')

	#print(transcripts.get('NM_000426.3'))

	egl = egl_variants(filename = 'data/EmVClass.2020-Q3.csv')
	#egl.to_csv('test.tsv',index=False,sep='\t')

	merge_datasets(clinvar,lovd,egl,output='LAMA2_pathogenic_variants.tsv')
	#merge_datasets(output = './test.vcf')
