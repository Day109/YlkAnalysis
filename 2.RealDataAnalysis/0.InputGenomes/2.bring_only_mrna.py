import glob, os

os.system('rm *_mRNA.gtf')

flist = glob.glob('*.gtf')

for f in flist:
	print(f)
	infile = open(f,"r")
	outfile = open(f.replace('.gtf','_mRNA.gtf'),"w")
	for sLine in infile:
	
		if not sLine.startswith("#"):
			sList = sLine.strip().split("\t")
			sAttribute = sList[8]
		
			if 'gene_biotype "protein_coding"' in sAttribute:
				outfile.write(sLine)

	infile.close()
	outfile.close()
			
