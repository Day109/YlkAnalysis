# YlkAnalysis
This page contains the custom python scripts used in an article "Detecting 
interspecies selection signal using population genomes and sequence 
divergences from the ancestor".


The programs used in 1.Simulation analysis includes:

- FastMl(v3.11)
- vcf2fasta
- SFS_CODE
- Prank(v170427)

and 2.Real data analysis includes:

- FastMl(v3.11)
- Prank(v170427)

These programs need to be in "1.SimulationAnalysis/Program" directory and
"2.RealDataAnalysis/Program" directory, respectively.

Moreover, to run the python script, the following modules need to be installed:

- pandas
- scipy
- Bio.Seq


## Running the pipeline ###

### Simulation data analysis

1.Simulation analysis consists of simplified population history of Great Ape lineage. 
To run this analysis, first move to "1.SimulationAnalysis/TestScenario" and type:

python 1.2.SfsCodeSimulation.py


### Real data analysis

2.Real data analysis is performed using the Great Ape genome data available 
from Ensembl (Whole genome fasta, gtf and vcf files). To run the pipeline, 
move to "2.RealDataAnalysis/" and unzip real data by typing:

cat YlkRealData.tar.gz.parta* | tar zxv;

Then, go to "2.RealDataAnalysis/Scripts/" and run the python scripts 
"1.Preprocessing.py", "2.AncestralMK.py" 
by typing:

python 1.Preprocessing.py # Don't have to run this line since the processed 
                                                        files are already in the main directory;

# Example commandline for 2.AncestralMK.py can be found in "ExampleCommandline.txt"
