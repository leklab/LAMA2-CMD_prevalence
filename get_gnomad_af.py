import requests

def get_gene_details_b38(gene_name):

  query = """
  {
    gene(gene_symbol: "%s", reference_genome: GRCh38){
      variants(dataset: gnomad_r3){
        chrom
        pos
        ref
        alt
        consequence_in_canonical_transcript
        consequence
        hgvsp
        hgvs
        genome{
          ac
          an
          populations{
            id
            ac
            an
          }
        }
      }
    }
  }""" % (gene_name)

  #print(query)

  res = requests.post('http://gnomad.broadinstitute.org/api', json={'query': query})

  if res.ok:
    #print res.json()
    return res.json()

  else:
    return({'data': 'error'}) 


def get_gene_details_b37(gene_name):

  query = """
  {
    gene(gene_symbol: "%s", reference_genome: GRCh37){
      variants(dataset: gnomad_r2_1){
        chrom
        pos
        ref
        alt
        consequence_in_canonical_transcript
        consequence
        hgvsp
        hgvs
        genome{
          ac
          an
          populations{
            id
            ac
            an
          }
        }
        exome{
          ac
          an
          populations{
            id
            ac
            an
          }
        }
      }
    }
  }""" % (gene_name)

  #print(query)

  res = requests.post('http://gnomad.broadinstitute.org/api', json={'query': query})

  if res.ok:
    #print res.json()
    return res.json()

  else:
    return({'data': 'error'}) 



if __name__ == '__main__':

  '''
  Gnomad R2_1 in B37 coordinates
  Contains exomes and genomes
  '''

  with open('LAMA2_hg19_hg38_liftover.txt') as f:
    lines = [line.rstrip() for line in f]

  #print(lines[2])

  data = get_gene_details_b37('LAMA2')


  '''
  for x in data['data']['gene']['variants']:
    if x['genome'] and x['exome']:
      print("%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d" % (x['chrom'],x['pos'],x['ref'],x['alt'],x['consequence'],x['hgvs'],x['genome']['ac'],x['genome']['an'],x['exome']['ac'],x['exome']['an']))

    elif x['genome']:
      print("%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\tNA\tNA" % (x['chrom'],x['pos'],x['ref'],x['alt'],x['consequence'],x['hgvs'],x['genome']['ac'],x['genome']['an']))

    elif x['exome']:
      print("%s\t%d\t%s\t%s\t%s\t%s\tNA\tNA\t%d\t%d" % (x['chrom'],x['pos'],x['ref'],x['alt'],x['consequence'],x['hgvs'],x['exome']['ac'],x['exome']['an']))
  '''

  #print("Liftover file length: %d\ngnomAD data length: %d" %(len(lines),len(data['data']['gene']['variants'])))


  pop_af_order = "ALL_AC\tALL_AN\tAFR_AC\tAFR_AN\tAMR_AC\tAMR_AN\tASJ_AC\tASJ_AN\tEAS_AC\tEAS_AN\tFIN_AC\tFIN_AN\tNFE_AC\tNFE_AN\tOTH_AC\tOTH_AN\tSAS_AC\tSAS_AN"
  #print("CHROM\t")
  
  gnomad_af_lookup = {}

  for i in range(len(data['data']['gene']['variants'])):
    x = data['data']['gene']['variants'][i]
    
    hg38_pos = int(lines[i].split('-')[1])
    variant_key = "%s-%d-%s-%s" % (x['chrom'],hg38_pos,x['ref'],x['alt'])

    if x['genome'] and x['exome']:
      genome_exome_pop_af = "Exome_r2_1,Genome_r2_1\t%s\t%s\t%d\t%d" % (x['consequence'],x['hgvs'],x['genome']['ac']+x['exome']['ac'], x['genome']['an']+x['exome']['an'])

      for j in range(len(x['genome']['populations'])):
        genome_pop = x['genome']['populations'][j]
        exome_pop = x['exome']['populations'][j]

        genome_exome_pop_af += "\t%d\t%d" % (genome_pop['ac']+exome_pop['ac'],genome_pop['an']+exome_pop['an'])


      gnomad_af_lookup[variant_key] = genome_exome_pop_af


    elif x['genome']:

      genome_pop_af = "Genome_r2_1\t%s\t%s\t%d\t%d" % (x['consequence'],x['hgvs'],x['genome']['ac'], x['genome']['an'])

      for j in range(len(x['genome']['populations'])):
        genome_pop = x['genome']['populations'][j]
        genome_pop_af += "\t%d\t%d" % (genome_pop['ac'],genome_pop['an'])


      gnomad_af_lookup[variant_key] = genome_pop_af



    elif x['exome']:

      exome_pop_af = "Exome_r2_1\t%s\t%s\t%d\t%d" % (x['consequence'],x['hgvs'],x['exome']['ac'], x['exome']['an'])

      for j in range(len(x['exome']['populations'])):
        exome_pop = x['exome']['populations'][j]
        exome_pop_af += "\t%d\t%d" % (exome_pop['ac'],exome_pop['an'])


      gnomad_af_lookup[variant_key] = exome_pop_af


'''
Gnomad R3 in B38 coordinates
Gnomad R3 is genomes only
'''


data = get_gene_details_b38('LAMA2')
gnomad_variant_keys = gnomad_af_lookup.keys()

for x in data['data']['gene']['variants']:

  
  variant_key = "%s-%d-%s-%s" % (x['chrom'],x['pos'],x['ref'],x['alt'])
  
  if variant_key in gnomad_variant_keys:
    continue

  genome_pop_af = "Genome_r3\t%s\t%s\t%d\t%d" % (x['consequence'],x['hgvs'],x['genome']['ac'], x['genome']['an'])

  pop_index_order = [8, 1, 7, 4, 2, 9, 0, 6]

  #for j in range(len(x['genome']['populations'])):
  for j in pop_index_order:

    genome_pop = x['genome']['populations'][j]

    #if genome_pop['id'] == 'ami' or genome_pop['id'] == 'mid':
      #continue

    genome_pop_af += "\t%d\t%d" % (genome_pop['ac'],genome_pop['an'])


  #print("%s\t%d\t%s\t%s\t%s" % (x['chrom'],x['pos'],x['ref'],x['alt'],genome_pop_af))
  gnomad_af_lookup[variant_key] = genome_pop_af



#6 129478729 clinvar_558356  C CAAAG . . Pathogenic|dup



pop_af_order = "ALL_AC\tALL_AN\tAFR_AC\tAFR_AN\tAMR_AC\tAMR_AN\tASJ_AC\tASJ_AN\tEAS_AC\tEAS_AN\tFIN_AC\tFIN_AN\tNFE_AC\tNFE_AN\tOTH_AC\tOTH_AN\tSAS_AC\tSAS_AN"
vcf_header_str = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
#print("CHROM\t")


print("%s\tgnomAD_Source\tAnnotation\thgvs\t%s" %(vcf_header_str,pop_af_order))


lama2_patho_variants = open('LAMA2_pathogenic_variants.tsv')
gnomad_variant_keys = gnomad_af_lookup.keys()
gnomad_known_pathogenic = {}

for line in lama2_patho_variants:
  fields = line.strip('\n').split('\t')

  if len(fields[3]) > 10 or len(fields[4]) > 10:
    continue

  variant_key = "%s-%d-%s-%s" % (fields[0],int(fields[1]),fields[3],fields[4])

  if variant_key in gnomad_variant_keys:
    print("%s\t%s" %(line.strip(),gnomad_af_lookup[variant_key]))
    gnomad_known_pathogenic[variant_key] = 1




'''
Grab novel loss of function variants

'''


gnomad_known_pathogenic_keys = gnomad_known_pathogenic.keys()


lof_variants = ['frameshift_variant','splice_acceptor_variant','splice_donor_variant','stop_gained']


for k in gnomad_variant_keys:
  fields = gnomad_af_lookup[k].split('\t')
  annotation = fields[1]
  #print(annotation)
  key_fields = k.split('-')

  if (annotation in lof_variants) and (not k in gnomad_known_pathogenic_keys):
  #if (annotation in lof_variants):
    vcf_line = "%s\t%d\t.\t%s\t%s\t.\t.\tNovel_gnomAD_LoF" % (key_fields[0],int(key_fields[1]),key_fields[2],key_fields[3])
    print("%s\t%s" %(vcf_line,gnomad_af_lookup[k]))












