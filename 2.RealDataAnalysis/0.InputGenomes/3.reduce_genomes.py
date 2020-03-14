import os, glob

for GtfFile in sorted(glob.glob('*_mRNA.gtf')):
	print('Parsing Chromosome list in '+GtfFile)
	ChromList = [];Species = GtfFile.split('.')[0]
	fin = open(GtfFile,'r')
	for line in fin:

		Chrom = line.strip().split('\t')[0]
		if not Chrom in ChromList:ChromList.append(Chrom)

	fin.close()

	FastaFile = glob.glob(Species+'*.dna.toplevel.fa')[0]
	print('Reducing '+FastaFile)
	fout = open(FastaFile.replace('.dna.toplevel.fa','.Reduced.fa'),'w')
	fin = open(FastaFile,'r');Flag = 0
	for line in fin:

		if '>' in line[0]:
			Chrom = line.strip().split()[0].replace('>','')
			if Chrom in ChromList:
				Flag = 1
				fout.write(line)
			else:
				Flag = 0
		else:
			if Flag == 1:fout.write(line)

	fin.close();fout.close()

