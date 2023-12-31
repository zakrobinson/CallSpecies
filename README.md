# CallSpecies.py
A python script for species indentification based on multilocus genotypes of species-informative markers.\
<b>Authors:</b> Zak Robinson and Jeff Stephenson; Robinson et al. 2023 <i>In Review</i> \
<b>Contact:</b> zrobinson@critfc.org

This script was written to fit seemlessly into existing GTseq genotyping pipelines \(e.g., <https://github.com/zakrobinson/GTseq-pipeline>\).
Please feel free to reach out if you would like assistance modifying these to better fit your work flow. 

## Quick Start: Example Run 

### Using Real Data
```
./Callspecies.py --help 

./CallSpecies.py --SpeciesSeq SpeciesSeq.csv --inGENO TESTDATA-GenosNF.csv --outGENO TESTDATA-GenosNF_CALLED.csv 
```
### Using Simulated Genotypic Data
```
./MakeMissingGenos.py --SpeciesSeq SpeciesSeq.csv --iters 5 --outGENO missingTEST-Genos.csv

./CallSpecies.py --SpeciesSeq SpeciesSeq.csv --inGENO missingTEST-Genos.csv --outGENO missingTEST-Genos_CALLED.csv 
```
Refer to MANUAL_CALLSPECIES.pdf for detailed description

