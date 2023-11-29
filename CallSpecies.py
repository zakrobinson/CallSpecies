#! /usr/bin/env python3

# CallSpecies.py Version 1.0, March 2023, Zak Robinson and Jeff Stephenson
# Contact: zrobinson@critfc.org




import argparse 
import copy
des="Add species calls to GenoCompile.pl derived genotype file\n"
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Add species calls to GenoCompile_v5.1.py derived genotype file",usage="%(prog)s --SpeciesSeq --inGENO --outGENO [--buffMP] [--pruneMP] [--colSTRT] [--outScoreMat] [--outRepGeno] [-h --help]",add_help=True)
    parser.add_argument('--SpeciesSeq', dest="species_seq", metavar='',required=True, help="SpeciesSeq file path: Locus,A1,A2,SpeciesA1,SpeciesA2,Weight")
    parser.add_argument('--inGENO', dest="inGENO", metavar='',required=True, help="File path for wide-format, comma-separated genotypic input file i.e., GenoCompile.pl derived file")
    parser.add_argument('--outGENO', dest="outGENO", metavar='',required=True, help="File path for genotypic output file with species calls fields added")
    parser.add_argument('--thresHET',dest="thresHET",metavar='',required=False,default=0.06,type=float,help="Threshold level of individual heterozygosity at species markers after which the species call is flagged with '"'ReviewHET;'"' default: %(default)s")
    parser.add_argument('--thresMS',dest="thresMS",metavar='',required=False,default=0.5,type=float,help="Threshold proportion of maximum absolute score that atleast one species must obtain to make a species call for a sample. default: %(default)s")
    parser.add_argument('--buffMP',dest="buffMP",metavar='',required=False,default=0.001,type=float,help="The buffer range around the highest proportion of maximum possible score (given missing data and heterozygotes) to consider a species as a candidate for a sample. default: %(default)s")
    parser.add_argument('--pruneMS',dest="pruneMS",metavar='',required=False,default=0.34,type=float,help="Threshold for proportion of maximum absolute score to retain a species as a candidate for a sample.  default: %(default)s")
    parser.add_argument('--colSTRT',dest="colSTRT",metavar='',required=False,default=5,type=int,help="1-based column position to insert two species call fields in --outGENO.  default: %(default)s")        
    parser.add_argument('--outScoreMat',dest="outScoreMat",metavar='',required=False, help="If specified, output file path for species score matrix (Optional)")
    parser.add_argument('--outRepGeno',dest="outRepGeno",metavar='',required=False, help="If specified, output file path for representative genotype per species (Optional)")
    args = parser.parse_args()
    if args.thresHET < 0 or args.thresHET > 1:
      parser.error('--thresHET must be between 0 and 1\n')
    if args.thresMS < 0 or args.thresMS > 1:
      parser.error('--thresMS must be between 0 and 1\n')
    if args.buffMP < 0 or args.buffMP > 1:
      parser.error('--buffMP must be between 0 and 1\n')
    if args.pruneMS < 0 or args.pruneMS > 1:
      parser.error('--pruneMS must be between 0 and 1\n')
    
      

######## Function to wrap up the scoring of a single individual from Prog90 ##############
def scoreInd(ind_Gdict):
  global SpeciesScore_Dict
  global SpeciesSeq_dict
  ind_score=copy.deepcopy(SpeciesScore_Dict) # create a copy the species score dictionary for scoring a single individual
  ind_hetCount=0 # initiate a count of heterozygotes at zero. Reported as a coarse index of potential hybridity or contamination. 
  nloc = 0
  for loc in ind_Gdict.keys():
    geno=ind_Gdict[loc].split(":")
    if ind_Gdict[loc] in ["0:0","0"]:
      continue # skipping missing data for now. 
    elif geno[0] != geno[1]:
      nloc+=1
      ind_hetCount+=1 # Not scoring hets; just counting them.
    elif geno[0] == SpeciesSeq_dict[loc][0]: # if the first allele matches the first allele in SpeciesSeq add weight for SpeciesA1 and substract for SpeciesA2.
      nloc+=1
      for s in SpeciesSeq_dict[loc][2].split(";"):
        ind_score[s]+=int(SpeciesSeq_dict[loc][4])
      for s in SpeciesSeq_dict[loc][3].split(";"):
        ind_score[s]-=int(SpeciesSeq_dict[loc][4])
    elif geno[0] == SpeciesSeq_dict[loc][1]: # if the first allele matches the second allele in SpeciesSeq add weight for SpeciesA2 and substract for SpeciesA1.
      nloc+=1
      for s in SpeciesSeq_dict[loc][2].split(";"):
        ind_score[s]-=int(SpeciesSeq_dict[loc][4])
      for s in SpeciesSeq_dict[loc][3].split(";"):
        ind_score[s]+=int(SpeciesSeq_dict[loc][4])
    else:
      pass # maybe add some error checking here or above. Anticipate key errors from Species Seq?
  hetText=str(ind_hetCount) + ";" + str(nloc) 
  hetfloat= 0 if nloc==0 else ind_hetCount/nloc
  return([ind_score,hetText,hetfloat])
##########################################################################################
#
#
#### Function to calculate max possible score given individual missingness and het calls #
def max_score_missing(ind_Gdict):
  global SpeciesScore_Dict
  global SpeciesSeq_dict
  mSMiss_dict=copy.deepcopy(SpeciesScore_Dict) # Create an empty score card
  for loc in SpeciesSeq_dict.keys():
    if loc not in ind_Gdict.keys():
      continue
    loc_specs=SpeciesSeq_dict[loc][2].split(";")
    loc_specs.extend(SpeciesSeq_dict[loc][3].split(";")) # Get a list of all species corresponding to this locus
    geno=ind_Gdict[loc].split(":")
    if ind_Gdict[loc] not in ["0:0","0"] and geno[0]==geno[1]: # if the individual is not missing data or heterozygous at this locus proceed with scoring
      for s in loc_specs:
        mSMiss_dict[s]+=int(SpeciesSeq_dict[loc][4])
  return(mSMiss_dict)
#########################################################################################
#
#
#### Read in SpeciesSeq file as dictionary {Locus:[A1,A2,SpeciesA1,SpeciesA2,Weight]} ####
SpeciesSeq_dict={}
with open(args.species_seq,"r") as SpeciesSeq:
  next(SpeciesSeq) #skip header
  for line in SpeciesSeq:
    info=line.strip().split(",")
    SpeciesSeq_dict[info[0]]=info[1:]
##########################################################################################
#
#
######## Unique Species from SpeciesSeq file. Create generic "score card" dict ###########
speciesConsidered=[]
for loc in SpeciesSeq_dict.keys():
  loc_specs=SpeciesSeq_dict[loc][2].split(";")
  loc_specs.extend(SpeciesSeq_dict[loc][3].split(";"))
  [speciesConsidered.append(i) for i in loc_specs if i not in speciesConsidered]
#generic score card below
SpeciesScore_Dict={spec:0 for spec in speciesConsidered} # generic species scoring dictionary initiated at score of zero. 
##########################################################################################
#
#
############# Calculate Absolute Maximum Possible Score per Species #######################
# Max score defined by SpeciesSeq file
absolute_maxS_dict=copy.deepcopy(SpeciesScore_Dict)
for loc in SpeciesSeq_dict.keys():
  loc_specs=SpeciesSeq_dict[loc][2].split(";")
  loc_specs.extend(SpeciesSeq_dict[loc][3].split(";"))
  for s in loc_specs:
    absolute_maxS_dict[s]+=int(SpeciesSeq_dict[loc][4])
##########################################################################################
#
#
########### Create representative genotypes according to SpeciesSeq ######################

if args.outRepGeno:
  specREPout=open(args.outRepGeno,"w")
  specREPout.write("Species,"+",".join(SpeciesSeq_dict.keys())+"\n")
perfectSpecies_dict={}  
for s in speciesConsidered:
  perfectGeno_dict={}
  for loc in SpeciesSeq_dict.keys():
    if s in SpeciesSeq_dict[loc][2].split(";"):
      geno=":".join([SpeciesSeq_dict[loc][0]]*2) # Generate a homozygous genotype
      perfectGeno_dict[loc]=geno
    elif s in SpeciesSeq_dict[loc][3].split(";"):
      geno=":".join([SpeciesSeq_dict[loc][1]]*2)
      perfectGeno_dict[loc]=geno
    else:
      perfectGeno_dict[loc]="0:0"
  if args.outRepGeno:
    specREPout.write(s + "," + ",".join(perfectGeno_dict.values())+"\n")
  perscore,junkhetT,junkhetF=scoreInd(ind_Gdict=perfectGeno_dict)
  perMAXscore={k:perscore[k]/absolute_maxS_dict[k] for k in perscore.keys()}
  maxscore=max(perMAXscore.values()) # maximum score among all species for individual line
  species_call=[k for k in perMAXscore.keys() if perMAXscore[k]>(maxscore-0.001)] # all canididates with max score are reported in list. Note the offset to avoid issues with float comparisons
  species_call=";".join(species_call) # report all candidate species as colon delimited string, hopefully one species is called
  perfectSpecies_dict[s]=",".join([str(perscore[k])+";"+str(absolute_maxS_dict[k]) for k in perscore.keys()])+","+str(junkhetT)+","+species_call+"\n"

if args.outRepGeno:
  specREPout.close()

##########################################################################################
#
#
########## Open outScoreMat connection and write headerline and max score line ###########
if args.outScoreMat:  # only write score matrix if file is specified
  # Prep species score matrix header line. Open file connection.
  Species_Score_HL=list(SpeciesScore_Dict.keys()) # Determine the species we are scoring from SpeciesScore_Dict at top of script.
  #add fields
  Species_Score_HL.insert(0,"Name")
  Species_Score_HL.insert(0,"Type")
  Species_Score_HL.extend(["SpeciesHET","SpeciesCall"])
  # open file connection and write headerline. 
  specMAT_out=open(args.outScoreMat,"w")
  specMAT_out.write(",".join(Species_Score_HL)+"\n")
  specMAT_out.write("SUMMARY,MAX_VALUES"+","+",".join([str(i) for i in absolute_maxS_dict.values()])+",NA,NA\n")
  for k in perfectSpecies_dict.keys():
    specMAT_out.write("SPECIES_REP," + k +"," + perfectSpecies_dict[k])
##########################################################################################
#
#
########################Open Input amd Output ProgGenos90 ################################
Prog90=open(args.inGENO,"r")
header=Prog90.readline() # read Prog90 header line 
header_list=header.strip().split(",") # Parse CSV header
orig_header=header_list.copy() # store originial header for reference in scoring function.
locus_index=[i for i in range(len(header_list)) if header_list[i] in SpeciesSeq_dict.keys()] # Where those Species Markers At. Will need to reference
## Print Info about Loci  in SpeciesSeq and ProgGenos90 input
nloci_specSeq=len(SpeciesSeq_dict.keys()) 
nloci_found=len(locus_index)
markers_notfound = [i for i in SpeciesSeq_dict.keys() if i not in header_list]
print(f'\nSpeciesSeq file contains : {nloci_specSeq} loci\nInput GenoCompile file contains: {nloci_found} loci in common\nLoci not shared:\n\n{markers_notfound}\n') # Just starting with this. Could certainly be more informative.
## Update header line with new fields
header_list.insert(args.colSTRT-1,"SpeciesCall")
header_list.insert(args.colSTRT,"SpeciesHET")
## Open new prog90 output w/ species call and write header line. 
NewProg90_out=open(args.outGENO,"w")
NewProg90_out.write(",".join(header_list)+"\n")
##########################################################################################
#
#
###################### Actually make the calls ###########################################
for line in Prog90:     # header line is already read in
  ind_list=line.strip().split(",") # get an individual line 
  ind_geno_dict={orig_header[i]:ind_list[i] for i in locus_index} # This creates an individual genotype dictionary {Locus:Genotype} just considering those shared between SpeciesSeq and GenoCompile input
  ind_score,het_text,het_float=scoreInd(ind_Gdict=ind_geno_dict) # Calculate individual species score dictionary, formatted het count "hets:non-missing genotypes", and heterozygosity of loci shared between SpeciesSeq and GenoCompile input
  genosucc=int(het_text.split(";")[1])/nloci_found if nloci_found>0 else 0
  indMAXdict=max_score_missing(ind_Gdict=ind_geno_dict) # Maximum score per species given missing data and het calls in this individual
  indMAXdict={k:(v if v!=0 else -999) for (k,v) in indMAXdict.items()} # remove the possibility of division by zero below.
  
  ind_perABSMAX={k:ind_score[k]/absolute_maxS_dict[k] for k in ind_score.keys()} # Calculate proportion of absolute max score for each individual.
  ind_perRELMAX={k:ind_score[k]/indMAXdict[k] for k in ind_score.keys()} # Calculate proportion of possible max score for the individual considering missing data and het calls.
  maxscoreABS=max(ind_perABSMAX.values()) # maximum absolute score among all species for individual line
  #maxscoreREL=max(ind_perRELMAX.values()) # maximum relative score among all species for an individual line
  
  if maxscoreABS < args.thresMS: # Despite the straightforwardness of missing data this allows more flexibility. Not all loci genotype for each species etc. It corresponds to missing data for each species when weights are 1.
    species_call="NoCall"
    if het_float > args.thresHET and genosucc >= 0.5: # This is the only time genotype success enters the picture directly. We dont want to flag review hets for failed samples. That fail the initial maxscore filter.
      species_call="ReviewHET;"+species_call
  else:
    maxscoreREL=max([v for (k,v) in ind_perRELMAX.items() if ind_perABSMAX[k]>=args.thresMS]) # maximum relative score among all species for an individual line (given missing data and hets); given the --thresMS is achieved.
    species_call=[k for k in ind_perRELMAX.keys() if ind_perRELMAX[k]>=(maxscoreREL-args.buffMP)] # all canididates with max score are reported in list. Note the offset, it avoid issues with float comparisons and permits approximate ties.
    if args.pruneMS:
      species_call=[spec for spec in species_call if ind_perABSMAX[spec]>=args.pruneMS] # The species called must achieve the --pruneMP max absolute score.
      species_call=["NoCall"] if not species_call else species_call # If none pass --thresMS, its a NoCall

    species_call=";".join(species_call) # report all candidate species as colon delimited string, hopefully one species is called
    if het_float>args.thresHET:
      species_call="ReviewHET;"+ species_call

  if args.outScoreMat: # if outScoreMat argument is specified write to file
    specMAT_out.write("sample,"+ ind_list[0]+","+",".join([str(ind_score[k])+";"+ str(indMAXdict[k] if indMAXdict[k] >0 else 0) for k in ind_score.keys()])+","+ het_text +","+species_call+"\n") # write the species score matrix line.
  
  ind_list.insert(args.colSTRT-1,species_call) # add fields to prog90 individual line
  ind_list.insert(args.colSTRT,het_text)
  NewProg90_out.write(",".join(ind_list)+"\n") # Write line to new prog90 output. We didn't mess with original data other than using .insert - should be good. We even retain the terrible last comma
##########################################################################################
#
#
##################### Close Open File Connections#########################################
NewProg90_out.close()
Prog90.close()
if args.outScoreMat:
  specMAT_out.close()
##########################################################################################

print("Done Calling Species\n")
