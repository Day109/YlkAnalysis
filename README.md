This page contains the custom python scripts used in an article "Detecting inter-species selection signal using population genomes and sequence divergence from the ancestor".

The programs used in 1.Simulation analysis includes:

- FastMl(v3.11)
- vcf2fasta
- SLiM(v3.3)
- Prank(v170427)

and 2.Real data analysis includes:

- FastMl(v3.11)
- Prank(v170427)

To run the python script, the following modules need to be installed:

- pandas
- scipy
- Bio.Seq

Before you could run these two analyses, you need to install the same versions of programs listed above in "Program" directories of "1.SimulationAnalysis" and "2.RealDataAnalysis". Also, to run the entire analysis, you need to have 700 Gb of space.

1.Simulation analysis consists of three different simulations resembling the simplified population history of Great Ape lineage. To run this analysis, move to the directories: "MutRate_Low", "MutRate_Med" and "MutRate_High" and run the python scripts "1.SimulationPipeline.py" and "2.SummarizeResults.py" by simply typing:

python 1.SimulationPipeline.py
python 2.SummarizeResults.py

2.Real data analysis is performed using the Great Ape genome data available from Ensembl (Whole genome fasta, gtf and vcf files). To run the pipeline, go to "Scripts" directory in "2.RealDataAnalysis" and run the python scripts "1.Preprocessing.py", "2.1.MkYlkPipeline.py" and "2.2.MkYlk_FastMl.py" by typing:

python 1.Preprocessing.py
python 2.1.MkYlkPipeline.py
python 2.2.MkYlk_FastMl.py

