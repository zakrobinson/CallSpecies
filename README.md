# CallSpecies.py
A python script for species indentification based on multilocus genotypes of species-informative markers.\
<b>Author:</b> Zak Robinson; Robinson et al. 2023 <i>In Review</i> \
<b>Contact:</b> zrobinson@critfc.org 

## Quick Start: Example Run 

### Using Real Data

./Callspecies.py --help 

./CallSpecies.py --SpeciesSeq SpeciesSeq.csv --inGENO TESTDATA-GenosNF.csv --outGENO TESTDATA-GenosNF_CALLED.csv 

### Using Simulated Genotypic Data

./MakeMissingGenos.py --SpeciesSeq SpeciesSeq.csv --iters 5 --outGENO missingTEST-Genos.csv

./CallSpecies.py --SpeciesSeq SpeciesSeq.csv --inGENO missingTEST-Genos.csv --outGENO missingTEST-Genos_CALLED.csv 

Refer to MANUAL_CALLSPECIES.pdf for detailed description

