#! /usr/bin/env python3

# MakeSharedPoly.py Version 1.0, Zak Robinson
# Contact: zrobinson@critfc.org

import argparse
import random
from scipy.stats import binom

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simulate Genotypes for Species Informative Markers with Shared Polymorphisms",usage="./%(prog)s --SpeciesSeq --outGENOS --iters --EffSpecies --nLociPoly --altFREQ [-h --help]",add_help=True)
    parser.add_argument('--SpeciesSeq', dest="species_seq", metavar='',required=True, help="SpeciesSeq file path: Locus,A1,A1,SpeciesA1,SpeciesA2,Weight")
    parser.add_argument('--outGENOS', dest="out_file", metavar='',required=True, help="File path for wide-format, comma-separated genotypic output file i.e., GenoCompile_v5.1.py formatted file")
    parser.add_argument('--iters', dest="iters", metavar='',required=False,default=100, type=int, help="The number of simulated individuals for each species, nLociPoly, and altFREQ. Values much greater than 1000 are likely unnecessary and create a large file default: %(default)s")
    parser.add_argument('--EffSpecies', dest="EffSpecies", metavar='',required=False,default='all', type=str, help="A species to be affected by the shared polymorphism or all to affect all species default: %(default)s")
    parser.add_argument('--nLociPoly', dest="nLociPoly", metavar='',required=False,default='1,2,5',type=str, help="The number of loci with a shared polymorphism among species. Integer or comma-separated string of integers default: %(default)s")
    parser.add_argument('--altFREQ', dest="altFREQ", metavar='',required=False, default='0.01,0.05,0.1', type=str, help="The frequency of the allele not associated with the species. Float value or comma-separated string of float values default: %(default)s")
    args = parser.parse_args()


try:
  altfreqs=args.altFREQ.split(",")
  altfreqs=[float(i) for i in altfreqs]
except ValueError:
  raise Exception("--altFREQ expects a float value or comma-seperated float values e.g., 0.01,0.05")

try:
  numpoly=args.nLociPoly.split(",")
  numpoly=[int(i) for i in numpoly]
except ValueError:
  raise Exception("--nLociPoly expects non-zero integer or comma-separated integer values e.g., 1,10")


if sum([i<1 for i in numpoly]) > 0:
  raise Exception("--nLociPoly expects non-zero integer or comma-separated integer values e.g., 1,10")



SpeciesSeq_dict={}
with open(args.species_seq,"r") as SpeciesSeq:
  next(SpeciesSeq) #skip header
  for line in SpeciesSeq:
    info=line.strip().split(",")
    SpeciesSeq_dict[info[0]]=info[1:]
    SpeciesSeq_dict[info[0]][2] = SpeciesSeq_dict[info[0]][2].split(";")
    SpeciesSeq_dict[info[0]][3] = SpeciesSeq_dict[info[0]][3].split(";")

consideredSpecies=[]    
for v in SpeciesSeq_dict.values():
  spectmp=v[2].copy()
  spectmp.extend(v[3])
  for s in spectmp:
    if s not in consideredSpecies:
      consideredSpecies.append(s)

if args.EffSpecies == "all":
  args.EffSpecies=consideredSpecies
elif args.EffSpecies not in consideredSpecies:
  raise Exception("--EffSpecies not in Species Seq File")
else:
  args.EffSpecies=[args.EffSpecies]

spec_allele_dict={}
spec_loc_inform={}
for spec in consideredSpecies:
  spec_allele_dict[spec]={}
  spec_loc_inform[spec]=[]
  for loc in SpeciesSeq_dict.keys():
    spec_allele_dict[spec][loc]=[]
    if spec not in SpeciesSeq_dict[loc][2] and spec not in SpeciesSeq_dict[loc][3]:
      spec_allele_dict[spec][loc].extend(['0','0'])
    elif spec in SpeciesSeq_dict[loc][2]:
      spec_allele_dict[spec][loc].append(SpeciesSeq_dict[loc][0])
      spec_allele_dict[spec][loc].append(SpeciesSeq_dict[loc][1])
      spec_loc_inform[spec].append(loc)
    elif spec in SpeciesSeq_dict[loc][3]:
      spec_allele_dict[spec][loc].append(SpeciesSeq_dict[loc][1])
      spec_allele_dict[spec][loc].append(SpeciesSeq_dict[loc][0])
      spec_loc_inform[spec].append(loc)


header_out=["Sample","MaxMarkInform","nPoly","LociWpoly","FreqOfPoly","TrueSpecies","rep"]
header_out.extend(list(SpeciesSeq_dict.keys()))
indFMT='%0'+ str(len(str(args.iters))) + "d"
with open(args.out_file,'w') as outgenos:
  outgenos.write(",".join(header_out)+"\n")
  for s in args.EffSpecies:
    for af in altfreqs:
      for np in numpoly:
        if len(spec_loc_inform[s]) >= np:
          nlp=np
        else:
          continue # or nlp=len(spec_loc_inform[s]); But can't imagine the use case to continue generating data the same level of missing data. 
        for i in range(args.iters):
          locWpoly=random.sample(spec_loc_inform[s],nlp)
          indID="_".join((s,'%02d'%np,str(af),indFMT%(i+1)))
          out_line=[indID,str(len(spec_loc_inform[s])),str(nlp),";".join(locWpoly),str(af),s,str(i+1)]
          for loc in SpeciesSeq_dict.keys():
            if loc not in locWpoly:
              geno=spec_allele_dict[s][loc][0]+":"+ spec_allele_dict[s][loc][0]
              out_line.append(geno)
            else:
              aindex=[int(binom.rvs(n=1,p=af,size=1)),int(binom.rvs(n=1,p=af,size=1))]
              aindex.sort()
              geno=spec_allele_dict[s][loc][aindex[0]]+":"+ spec_allele_dict[s][loc][aindex[1]]
              out_line.append(geno)
          outgenos.write(",".join(out_line)+"\n")
print('done')
