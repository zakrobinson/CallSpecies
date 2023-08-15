#! /usr/bin/env python3
# MakeMissingGenos.py Version 1.0, Zak Robinson
# Contact: zrobinson@critfc.org

import random
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create All-levels of Missing Data for Species-Informative Markers",usage="./%(prog)s --SpeciesSeq --outGENOS --iters [-h --help]",add_help=True)
    parser.add_argument('--SpeciesSeq', dest="species_seq", metavar='',required=True, help="SpeciesSeq file path: Locus,A1,A1,SpeciesA1,SpeciesA2,Weight")
    parser.add_argument('--outGENOS', dest="outGENOS", metavar='',required=True, help="File path for wide-format, comma-separated genotypic output file i.e., GenoCompile_v5.1.py formatted file")
    parser.add_argument('--iters', dest="iters", metavar='',required=True, type=int, help="The number of random samples at each level of missing genotypic data for each species. Values greater than 1000 are likely unnecessary and create a large file")
    args = parser.parse_args()


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




header_out=["Sample","MaxMarkInform","nMarkers","PercGT","TrueSpecies","rep"]
header_out.extend(list(SpeciesSeq_dict.keys()))
indFMT='%0'+ str(len(str(args.iters))) + "d"
with open(args.outGENOS,"w") as fakeGenoOut:
  fakeGenoOut.write(",".join(header_out)+"\n")
  for s in consideredSpecies:
    n_inform=len(spec_loc_inform[s])
    for m in range(n_inform,0,-1):
      for i in range(args.iters):
        indID="_".join((s,'%03d'%m ,indFMT%(i+1)))
        include=random.sample(spec_loc_inform[s],m)
        outline=[indID,str(n_inform),str(m),str(m/n_inform),s,str(i+1)]
        for loc in SpeciesSeq_dict.keys():
          if loc not in include:
            outline.append("0:0")    
          else:
            geno=spec_allele_dict[s][loc][0]+":"+spec_allele_dict[s][loc][0]
            outline.append(geno)
        fakeGenoOut.write(",".join(outline)+"\n")

print("Done")


      


  
  
  
 
