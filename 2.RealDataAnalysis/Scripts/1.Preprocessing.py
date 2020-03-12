import glob, os, multiprocessing

BeingCompared = 'homo_sapiens'
ProgramDir = '../Program/'
MainDir = '../'

MsaThread = 50

TreeFile = '../GreatApe.tre'

def Parse_homology_information(FileName):

	## Parse homology information from Ensembl biomart export
	## biomart export needs to have column order as follows:
	## query gene, target gene, query protein, target protein, homology info
	## Input: Ensembl biomart export file name
	## Output: Homology information

	HomologyInfoList = []

	fin = open(FileName,'r')
	fin.readline()

	for line in fin:
		if len(line.strip().split('\t')) == 5:
			qGeneId,tGeneId,qProtId,tProtId,Hom = line.strip().split('\t')
			HomologyInfoList.append(line.strip().split('\t'))

	fin.close()

	return HomologyInfoList
			

def Parse_genes(SpeciesName):

	## Extract protein ids from gtf file
	## input: species name of gtf file
	## output: protein ids

	GeneList = {}

	GtfFile = glob.glob(MainDir+'0.InputGenomes/'+SpeciesName+'*_mRNA.gtf')
	fin = open(GtfFile[0],'r')
	for line in fin:
		if not '#' in line[0]:

			GeneType = line.strip().split('\t')[2]
			if 'CDS' in GeneType:
				P_id = line.strip().split('\t')[8].split('protein_id "')[1].split('"')[0]
				if not P_id in GeneList:
					GeneList.setdefault(P_id)

					if len(GeneList)%10000 == 0:
						print(str(len(GeneList))+' genes parsed')

	if len(GeneList)%10000 != 0:
		print(str(len(GeneList))+' genes parsed')
	
	fin.close()

	return GeneList


def Summarize_orthologue(Flist):

	## Summarize ortholog information for each gene
	## Input: list of Ensembl biomart export files
	## Output: OrthologGeneList.txt

	fout = open('OrthologGeneList.txt','w')
	fout.write('Query')

	QueryCdsList = Parse_genes(BeingCompared)
	One2oneOrthDic = {}

	for Cds in QueryCdsList:One2oneOrthDic.setdefault(Cds,'')
	for i in range(len(Flist)):

		Species = Flist[i].split('/')[-1].split('_vs_')[1]
		HomInfo = Parse_homology_information(Flist[i])

		fout.write('\t'+Species)

		for HomGene in HomInfo:
			qGeneId,tGeneId,qProtId,tProtId,Hom = HomGene

			if Hom == 'ortholog_one2one' and qProtId in QueryCdsList:
				One2oneOrthDic[qProtId] += '\t'+tProtId

		for Cds in QueryCdsList:
			if len(One2oneOrthDic[Cds].split('\t')) < i+2:
				One2oneOrthDic[Cds] += '\tNA'

	fout.write('\n')
	for Cds in QueryCdsList:
		fout.write(Cds+One2oneOrthDic[Cds]+'\n')
	fout.close()

	print('Ortholog genes summarized')


def Complement(Seq):
	Complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', \
	'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R', 'W': 'W', 'S': 'S', 'V': 'B', 'B': 'V', 'H': 'D', 'D': 'H'}
	Bases = list(Seq)
	Bases = [Complement[Base] for Base in Bases]
	return ''.join(Bases)
			

def Parse_GTF_file(GtfFile):
	
	## Parse GTF file
	## Input: GTF file name
	## Output: dictionary of GTF lines for each chromosome

	GtfDic = {}

	fin = open(GtfFile,'r')
	for line in fin:
		if not '#' in line[0]:
			GeneType = line.strip().split('\t')[2]
			if 'CDS' == GeneType:
				GeneId = line.strip().split('gene_id "')[1].split('"')[0]
				TransId = line.strip().split('transcript_id "')[1].split('"')[0]
				ProtId = line.strip().split('protein_id "')[1].split('"')[0]
				Chrom = line.strip().split('\t')[0].replace('chr','')
				Strand = line.strip().split('\t')[6]
				C_St = line.strip().split('\t')[3]
				C_End = line.strip().split('\t')[4]
				
				GtfDic.setdefault(Chrom,[])
				GtfDic[Chrom].append([GeneId,TransId,ProtId,Strand,C_St,C_End])

	fin.close()
	
	return GtfDic


def Parse_Cds_from_chromosome_sequence(GtfString, ChromSeq):

	## Parse Cds from GTF and Chromosome sequence
	## ProtSeq is DNA sequence of each protein coding gene
	## Input: Gtf string, Chromosome sequence
	## Output: DNA sequence of each protein coding gene

	ProtSeq = {}

	for line in GtfString:
		GeneId, TransId, ProtId, Strand, C_St, C_End = line
		ProtSeq.setdefault(TransId+'\t'+ProtId+'\t'+Strand,'')
		if Strand == '-':
			ProtSeq[TransId+'\t'+ProtId+'\t'+Strand] += Complement(ChromSeq[int(C_St)-1:int(C_End)].upper())[::-1]
		else:
			ProtSeq[TransId+'\t'+ProtId+'\t'+Strand] += ChromSeq[int(C_St)-1:int(C_End)].upper()

	return ProtSeq


def Generate_CdsFasta_sequence(GtfFile, ReferenceFile):

	## Parse Cds from GTF and Fasta files
	## Input: GTF file, Fasta file
	## Output: '2.CdsSeqeunces/' directory containing Cds fasta files

	Species = GtfFile.split('/')[-1].split('.')[0]
	print('Extracting Cds of '+Species)

	Cnt = 0;non_Cnt = 0

	if not '2.CdsSeqeunces' in os.listdir('../'):
		os.system('mkdir '+MainDir+'2.CdsSeqeunces')
	if not Species in os.listdir('../2.CdsSeqeunces'):
		os.system('mkdir '+MainDir+'2.CdsSeqeunces/'+Species)

	GtfDic = Parse_GTF_file(GtfFile)

	ChromSeq = '';RefChrom = ''
	ref_fin = open(ReferenceFile,'r');FirstLine = 1

	for line in ref_fin:
		if '>' in line[0]:
			if FirstLine == 1:
				pass
			else:
				if RefChrom in GtfDic:
					ProtSeq = Parse_Cds_from_chromosome_sequence(GtfDic[RefChrom], ChromSeq)

					for Cds in ProtSeq:
						TransId = Cds.split('\t')[0]
						ProtId = Cds.split('\t')[1]
						Strand = Cds.split('\t')[2]
						if len(ProtSeq[Cds]) % 3 == 0:
							fout = open(MainDir+'2.CdsSeqeunces/'+Species+'/'+TransId+'_'+ProtId+'.fa','w')
							fout.write('>'+Species+'\n'+ProtSeq[Cds])
							fout.close()
							Cnt += 1
						else:
							non_Cnt += 1

			FirstLine = 0
			RefChrom = line.strip().split()[0].replace('>','').replace('chr','')
			ChromSeq = ''

		else:
			ChromSeq += line.strip()

	ref_fin.close()

	#### Parsing Last Chromosome ####
	if RefChrom in GtfDic:
		ProtSeq = Parse_Cds_from_chromosome_sequence(GtfDic[RefChrom], ChromSeq)

		for Cds in ProtSeq:
			TransId = Cds.split('\t')[0]
			ProtId = Cds.split('\t')[1]
			Strand = Cds.split('\t')[2]
			if len(ProtSeq[Cds]) % 3 == 0: 
				fout = open(MainDir+'2.CdsSeqeunces/'+Species+'/'+TransId+'_'+ProtId+'.fa','w')
				fout.write('>'+Species+'\n'+ProtSeq[Cds])
				fout.close()
				Cnt += 1
			else:
				non_Cnt += 1

	print('Cds of '+Species+': '+str(Cnt))
	print('Cds of '+Species+' not in multiple of 3: '+str(non_Cnt))


def Multi_sample_GTF(Flist):

	## Generate Cds Fasta file from GTF and Fasta Files

	GtfFile = glob.glob(MainDir+'0.InputGenomes/'+BeingCompared+'*_mRNA.gtf')[0]
	ReferenceFile = glob.glob(MainDir+'0.InputGenomes/'+BeingCompared+'*.Reduced.fa')[0]
	Generate_CdsFasta_sequence(GtfFile, ReferenceFile)

	for f in Flist:
		Species = f.split('/')[-1].split('_vs_')[1]
		GtfFile = glob.glob(MainDir+'0.InputGenomes/'+Species+'*_mRNA.gtf')[0]
		ReferenceFile = glob.glob(MainDir+'0.InputGenomes/'+Species+'*.Reduced.fa')[0]

		Generate_CdsFasta_sequence(GtfFile, ReferenceFile)


def Merge_Cds_files(FileName):

	## Merging one2one ortholog genes
	## Input: Fasta files in '2.CdsSeqeunces/', 'OrthologGeneList.txt'
	## Output: '3.CdsMerged/' directory with Merged Fasta files

	if not '3.CdsMerged' in os.listdir(MainDir):
		os.system('mkdir '+MainDir+'3.CdsMerged')

	fin = open(FileName, 'r')
	SpeciesList = fin.readline().strip().split('\t')
	SpeciesList[0] = BeingCompared

	for line in fin:
		GeneList = line.strip().split('\t')

		if not 'NA' in GeneList:
			fout_string = []
			
			for i in range(len(SpeciesList)):

				SampleString = ''
				CdsFasta = glob.glob(MainDir+'2.CdsSeqeunces/'+SpeciesList[i]+'/*'+GeneList[i]+'*')
				
				if len(CdsFasta) == 1:

					fin = open(CdsFasta[0],'r')
					for line in fin:
						SampleString += line
					fin.close()

					SampleString += '\n'
					fout_string.append(SampleString)

			if len(fout_string) == len(SpeciesList):
				fout = open(MainDir+'3.CdsMerged/'+GeneList[0]+'.fa','w')
				for string in fout_string:
					fout.write(string)
		
				fout.close()

	fin.close()


def Execute_MSA(f):
	
	## Input: Merged Fasta file
	## Output: Aligned Fasta file

	GeneId = f.split('/')[-1].split('.')[0]

	systemstr = ProgramDir+'prank-msa/src/prank'+\
	' -d='+f+\
	' -o='+MainDir+'4.CdsAligned/'+GeneId+\
	' -f=fasta'+\
	' -quiet'+\
	' -F'+\
	' -codon'

	print(systemstr);os.system(systemstr)


def multi_MSA():

	## Multiprocessing of MSA
	## Input: Merged Fasta files in '3.CdsMerged/'
	## Output: Aligned Fasta files in '4.CdsAligned/'

	if not '4.CdsAligned' in os.listdir(MainDir):
		os.system('mkdir '+MainDir+'4.CdsAligned')

	MergedFlist = glob.glob(MainDir+'3.CdsMerged/*.fa')
	ReducedFlist = []
	for f in MergedFlist:
		GeneId = f.split('/')[-1].split('.')[0]
		if not GeneId+'.best.fas' in os.listdir(MainDir+'4.CdsAligned'):
			ReducedFlist.append(f)

	p = multiprocessing.Pool(MsaThread)
	p.map(Execute_MSA, ReducedFlist)
	p.close()
	p.join()


def Execute_MSA_WithTree(f):

	## Input: Merged Fasta file
	## Output: Aligned Fasta file

	GeneId = f.split('/')[-1].split('.')[0]

	systemstr = ProgramDir+'prank-msa/src/prank'+\
	' -d='+f+\
	' -o='+MainDir+'4.2.CdsAligned_WithAsr/'+GeneId+\
	' -f=fasta'+\
	' -t='+TreeFile+\
	' -showanc'+\
	' -showevents'+\
	' -F'+\
	' -once'+\
	' -quiet'+\
	' -codon'

	print(systemstr);os.system(systemstr)


def multi_MSA_with_ASR():

	## Multiprocessing of PRANK MSA with tree and ancestral sequence reconstruction
	## Input: Merged Fasta files in '3.CdsMerged/'
	## Output: Aligned Fasta files in '4.2.CdsAligned_WithAsr/'
	
	if not '4.2.CdsAligned_WithAsr/' in os.listdir(MainDir):
		os.system('mkdir '+MainDir+'4.2.CdsAligned_WithAsr/')
	
	MergedFlist = glob.glob(MainDir+'3.CdsMerged/*.fa')
	ReducedFlist = []
	for f in MergedFlist:
		GeneId = f.split('/')[-1].split('.')[0]
		if not GeneId+'.fa.best.fas' in os.listdir(MainDir+'4.2.CdsAligned_WithAsr'):
			ReducedFlist.append(f)

	p = multiprocessing.Pool(MsaThread)
	p.map(Execute_MSA_WithTree, ReducedFlist)
	p.close()
	p.join()


def Reorder_MSA(MsaFlist):

	## Reorder species into the order in Tree file
	## Input: Tree file, Aligned fasta files
	## Output: Reordered fasta files

	fin = open(TreeFile,'r')
	Tree = fin.readline()
	fin.close()

	SpeciesOrder = []	

	SpeciesList = Tree.strip().split(',')
	for element in SpeciesList:
		Species = element.replace('(','').replace(')','').replace(';','').replace(' #1','')
		SpeciesOrder.append(Species)

	for f in MsaFlist:

		GeneId = f.split('/')[-1].split('.')[0]

		SequenceDic = {}
		fin = open(f,'r')
		for line in fin:
			if '>' in line:
				Species = line.strip().replace('>','')
				SequenceDic.setdefault(Species,'')
			else:
				SequenceDic[Species] += line.strip().replace('N','-')
		fin.close()

		fout = open(f.replace('.best.fas','_reordered.fa'),'w')
		for Species in SpeciesOrder:
			fout.write('>'+Species+'\n'+SequenceDic[Species]+'\n')
		fout.close()


if __name__ == "__main__":

	Flist = glob.glob(MainDir+'1.HomologyFromEnsembl/*')

	Summarize_orthologue(Flist)
	
	Multi_sample_GTF(Flist)
	Merge_Cds_files('OrthologGeneList.txt')

	## Msa only
	multi_MSA()
	MsaFlist = glob.glob(MainDir+'4.CdsAligned/*.fas')
	Reorder_MSA(MsaFlist)

	## Msa with ancestral sequence reconstruction
	multi_MSA_with_ASR()
