This page contains the custom python scripts used in an article "Detecting 
inter-species selection signal using population genomes and sequence 
divergence from the ancestor".

## Requirements before running the pipeline ###

The programs used in 1.Simulation analysis includes:

- FastMl(v3.11)
- vcf2fasta
- SLiM(v3.3)
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

To run the entire analysis, you need to have 700 Gb of space.


## Running the pipeline ###

### Simulation data analysis

1.Simulation analysis consists of three different simulations resembling the 
simplified population history of Great Ape lineage. To run this analysis, 
first move to "1.SimulationAnalysis/" and unzip simulation data by:

python ExtractInputFiles.py

move to the directories: "MutRate_Low", "MutRate_Med" and "MutRate_High" and 
run the python scripts "1.SimulationPipeline.py" and "2.SummarizeResults.py" 
by simply typing:

python 1.SimulationPipeline.py
python 2.SummarizeResults.py


### Real data analysis

2.Real data analysis is performed using the Great Ape genome data available 
from Ensembl (Whole genome fasta, gtf and vcf files). To run the pipeline, 
move to "2.RealDataAnalysis/" and unzip real data by typing:

tar -zxvf YlkRealData.tar.gz

Then, go to "2.RealDataAnalysis/Scripts/" and run the python scripts 
"1.Preprocessing.py", "2.1.MkYlkPipeline.py" and "2.2.MkYlk_FastMl.py" 
by typing:

python 1.Preprocessing.py # Don't have to run this line since the processed 
							files are already in the main directory
python 2.1.MkYlkPipeline.py
python 2.2.MkYlk_FastMl.py
