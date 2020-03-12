import glob, os
from Bio.Seq import Seq
import pandas as pd

MainDir = '../'
TreeFile = MainDir+'GreatApe.tre'

def Create_GTF_df(GtfFile):

	## Parsing GTF file
	## Input: GTF file name
	## Output: Dataframe of GTF file and dictionary of gene length

	GtfDf = {};GeneLenDic = {}
	Header = ['Gene_ID','Tran_ID','Prot_ID','Exon_ID','Chrom','Strand','Start','End']
	for colname in Header:GtfDf.setdefault(colname,[])

	Species = GtfFile.split('/')[-1].split('.')[0]
	fin = open(GtfFile,'r');cnt = 0

	for line in fin:
		if not '#' in line[0]:
			Gtype = line.strip().split('\t')[2]
			if 'CDS' == Gtype:
				cnt += 1
				if cnt%100000==0:print('Parsed %s GTF lines: %d'%(Species,cnt))
				G_S_ID = line.strip().split('gene_id "')[1].split('"')[0]
				T_ID = line.strip().split('transcript_id "')[1].split('"')[0]
				if 'protein_id "' in line:P_ID = line.strip().split('protein_id "')[1].split('"')[0]
				else:P_ID = str(cnt)
				if 'exon_number "' in line:E_ID = line.strip().split('exon_number "')[1].split('"')[0]
				else:E_ID = str(cnt)
				Chrom = line.strip().split('\t')[0].replace('chr','')
				Strand = line.strip().split('\t')[6]
				E_St = int(line.strip().split('\t')[3])
				E_End = int(line.strip().split('\t')[4])

				GtfDf[Header[0]].append(G_S_ID)
				GtfDf[Header[1]].append(T_ID)
				GtfDf[Header[2]].append(P_ID)
				GtfDf[Header[3]].append(E_ID)
				GtfDf[Header[4]].append(Chrom)
				GtfDf[Header[5]].append(Strand)
				GtfDf[Header[6]].append(E_St)
				GtfDf[Header[7]].append(E_End)

				GeneLenDic.setdefault(P_ID,0)
				GeneLenDic[P_ID] += E_End-E_St+1

	fin.close()
	if cnt%100000!=0:print('Parsed %s GTF lines: %d'%(Species,cnt))

	GtfDf = pd.DataFrame(GtfDf, columns = Header)

	return GtfDf, GeneLenDic


def Convert_vcf_position_to_MSA():

	## Convert chromosomal positions in VCF file into Multiple sequence alignment position
	## Input: GTF file of protein-coding genes and VCF file
	## Output: Dataframe of variant loci from VCF file

	GtfList = glob.glob(MainDir+'0.InputGenomes/*_mRNA.Gtf')

	for Gtf in GtfList:
		Species = Gtf.split('/')[-1].split('.')[0]
		GtfDf, GeneLenDic = Create_GTF_df(Gtf)

		VcfDf = {};prev_chrom = ''
		cumulative_len_dic = {}
		vcf_Header = ['Chrom','Pos','New_pos','Ref','Alt','Gene_ID','Trans_ID','Prot_ID','Strand']
		for ColName in vcf_Header:VcfDf.setdefault(ColName,[])
		
		cnt = 0
		vcf_file = MainDir+'0.2.VcfFiles/'+Species+'.vcf'
		fin = open(vcf_file,'r')
		for line in fin:
			cnt += 1
			if cnt%1000==0:print('Parsing %s VCF line: %d'%(Species,cnt))
			if not line[0] == '#' and 'TSA=SNV' == line.split(';')[1]:
				chrom,pos,v_id,ref,alt = line.strip().split('\t')[:5]
				if chrom != prev_chrom:
					reduced_Gtf = GtfDf[GtfDf['Chrom']==chrom].sort_values(by=['Chrom','Start','End']).reset_index()

				if len(reduced_Gtf.index) != 0:	
					start = reduced_Gtf.iloc[0]["Start"];end = reduced_Gtf.iloc[0]["End"]
					g_id = reduced_Gtf.iloc[0]['Gene_ID']
					t_id = reduced_Gtf.iloc[0]['Tran_ID']
					p_id = reduced_Gtf.iloc[0]['Prot_ID']
					cumulative_len_dic.setdefault(p_id,0)

					if int(pos) > end:
						while int(pos) > end and len(reduced_Gtf.index) != 0:
							cumulative_len_dic[p_id] += end-start+1
							reduced_Gtf = reduced_Gtf.drop([0]).reset_index(drop=True)
							if len(reduced_Gtf.index) != 0:
								g_id = reduced_Gtf.iloc[0]['Gene_ID']
								t_id = reduced_Gtf.iloc[0]['Tran_ID']
								p_id = reduced_Gtf.iloc[0]['Prot_ID']
								cumulative_len_dic.setdefault(p_id,0)
								start = reduced_Gtf.iloc[0]["Start"]
								end = reduced_Gtf.iloc[0]["End"]

					elif int(pos) <= end and int(pos) >= start:
						if reduced_Gtf.iloc[0]['Strand'] == '+':
							new_pos = cumulative_len_dic[p_id] + int(pos) - start + 1 
						else:
							new_pos = GeneLenDic[p_id] - (cumulative_len_dic[p_id] + int(pos) - start)

						VcfDf[vcf_Header[0]].append(chrom)
						VcfDf[vcf_Header[1]].append(int(pos))
						VcfDf[vcf_Header[2]].append(new_pos)
						VcfDf[vcf_Header[3]].append(ref)
						VcfDf[vcf_Header[4]].append(alt)
						VcfDf[vcf_Header[5]].append(g_id)
						VcfDf[vcf_Header[6]].append(t_id)
						VcfDf[vcf_Header[7]].append(p_id)
						VcfDf[vcf_Header[8]].append(reduced_Gtf.iloc[0]['Strand'])
						
				prev_chrom = chrom
		fin.close()

		print('Parsing %s VCF line: %d'%(Species,cnt))

		VcfDf = pd.DataFrame(VcfDf, columns = vcf_Header)
		VcfDf.to_csv(MainDir+'0.2.VcfFiles/'+Species+'_VcfDf.txt',sep='\t')


def Complement(seq):
	Complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', \
	'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R', 'W': 'W', 'S': 'S', 'V': 'B', 'B': 'V', 'H': 'D', 'D': 'H', ',':','}
	bases = list(seq)
	bases = [Complement[base] for base in bases]
	return ''.join(bases)


def Count_P_n_D(Sseq,Tseq,CoordinateList,SvcfDf,GeneId):

	## Count polymorphism and divergence from aligned protein-coding sequences
	## Input: 1.Two protein-coding sequences to count divergences & 2.Coodinate list 
	##        conserved in Great ape species (gapless) & 3.Variant dataframe
	## Output: Countings of dN, dS, pN, pS

	dN = 0;dS = 0;pN = 0;pS = 0

	for base,i in CoordinateList:

		cst = int(i/3)*3;cend = cst+3
		OriCod = Sseq[cst:cend]

		### Remove Polymophism in gap ###
		while len(SvcfDf.index) != 0 and base > SvcfDf.iloc[0]['New_pos']:
			SvcfDf = SvcfDf.drop([0]).reset_index(drop=True)

		### No Divergence ###
		if Sseq[i] == Tseq[i]:

			### Polymorphism ###
			if len(SvcfDf.index) != 0 and base == (SvcfDf.iloc[0]['New_pos']-1):

				if SvcfDf.iloc[0]['Strand'] == '+':
					vrnts = SvcfDf.iloc[0]['Alt'].split(',')
				else:
					vrnts = Complement(SvcfDf.iloc[0]['Alt']).split(',')

				for vrnt in vrnts:
					VarCod = Sseq[cst:i]+vrnt+Sseq[i+1:cend]
					if Seq(OriCod).translate() == Seq(VarCod).translate():
						pS += 1
					else:
						pN += 1

				SvcfDf = SvcfDf.drop([0]).reset_index(drop=True)

		### Divergence ###
		elif Sseq[i] != Tseq[i]:

			### Polymorphism ###
			if len(SvcfDf.index) != 0 and base == (SvcfDf.iloc[0]['New_pos']-1):

				if SvcfDf.iloc[0]['Strand'] == '+':
					vrnts = SvcfDf.iloc[0]['Alt'].split(',')
				else:
					vrnts = Complement(SvcfDf.iloc[0]['Alt']).split(',')

				for vrnt in vrnts:
					VarCod = Sseq[cst:i]+vrnt+Sseq[i+1:cend]
					if Seq(OriCod).translate() == Seq(VarCod).translate():
						pS += 1
					else:
						pN += 1

				SvcfDf = SvcfDf.drop([0]).reset_index(drop=True)

				if not Tseq[i] in vrnts:
					DivCod = Sseq[cst:i]+Tseq[i]+Sseq[i+1:cend]
					if Seq(OriCod).translate() == Seq(DivCod).translate():
						dS += 1
					else:
						dN += 1

			### No Polymorphism ###
			else:
				DivCod = Sseq[cst:i]+Tseq[i]+Sseq[i+1:cend]
				if Seq(OriCod).translate() == Seq(DivCod).translate():
					dS += 1
				else:
					dN += 1

	return dN,dS,pN,pS,len(Sseq)


def Parse_reconstructed_seq(Fname):

	## Parse sequences from Fasta file

	SeqDic = {}
	fin = open(Fname,'r')
	for line in fin:
		if '>' in line[0]:
			fsid = line.strip().split(' ')[0].replace('>','')
			SeqDic.setdefault(fsid,'')
		else:SeqDic[fsid] += line.strip()
	fin.close()

	return SeqDic


def ParseHomologyDic(Species, FileName):

	## Parse one2one ortholog genes of a species

	HomDic = {}
	fin = open(FileName,'r')
	SpeciesList = pd.Series(fin.readline().split('\t'))
	SearchInd = SpeciesList[SpeciesList == Species].index[0]
	for line in fin:
		Genes = line.strip().split('\t')
		if not Genes[SearchInd] == 'NA':
			HomDic[Genes[0]] = Genes[SearchInd]

	return HomDic


def ParseGblockCoordinates(SeqDic,Species):

	## Parse gapless coordinate list of Multiple sequence alignment

	base = 0;CoordinateList = []
	for i in range(len(SeqDic[Species])):
		if SeqDic[Species][i] != '-':base += 1
		GapFlag = 0
		for SeqId in sorted(SeqDic.keys()):
			if '-' == SeqDic[SeqId][i]:GapFlag = 1
		if GapFlag == 0:CoordinateList.append([base-1,i])

	return CoordinateList


def Perform_MkYlk_test(Species,Target):

	## Perform Mk or Ylk test on a species of interest and a target
	## Input: Name of species and target sequences from MSA Fasta file
	## Output: Dataframe summarizing Mk or Ylk test

	IdDic = {'homo_sapiens':'homo_sapiens', 'pongo_abelii':'pongo_abelii', 
			 'gorilla_gorilla':'gorilla_gorilla', 'pan_troglodytes':'pan_troglodytes', 
			 'pan_paniscus':'pan_paniscus', '#4#':'1', '#3#':'2', '#2#':'3', '#1#':'4'}

	if not Species == 'homo_sapiens':
		HomDic = ParseHomologyDic(Species, 'OrthologGeneList.txt')

	VcfDf = pd.read_csv(MainDir+'0.2.VcfFiles/'+Species+'_VcfDf.txt',sep='\t',index_col=0)
	MkYlkDf = {};Header = ['Prot_ID','dN','dS','pN','pS','MsaLen']
	for colname in Header:MkYlkDf.setdefault(colname,[])

	AncMsa_list = sorted(glob.glob(MainDir+'4.2.CdsAligned_WithAsr/*.anc.fas'))
	for AncMsa in AncMsa_list:
		GeneId = AncMsa.split('/')[-1].split('.')[0]
		if not Species == 'homo_sapiens':GeneId = HomDic[GeneId]

		SeqDic = Parse_reconstructed_seq(AncMsa)
		CoodinateList = ParseGblockCoordinates(SeqDic,Species)
		reduced_VcfDf = VcfDf[VcfDf['Prot_ID']==GeneId]
		reduced_VcfDf = reduced_VcfDf.sort_values(by='New_pos').reset_index(drop=True)

		dN,dS,pN,pS,MsaLen = Count_P_n_D(SeqDic[Species],SeqDic[Target],CoodinateList,reduced_VcfDf,GeneId)

		MkYlkDf[Header[0]].append(GeneId)
		MkYlkDf[Header[1]].append(dN)
		MkYlkDf[Header[2]].append(dS)
		MkYlkDf[Header[3]].append(pN)
		MkYlkDf[Header[4]].append(pS)
		MkYlkDf[Header[5]].append(MsaLen)

	Species = IdDic[Species];Target = IdDic[Target]
	MkYlkDf = pd.DataFrame(MkYlkDf, columns=Header)
	MkYlkDf.to_csv(Species+'_'+Target+'_Prank-MkYlk.txt',sep='\t')


if __name__ == "__main__":

	YlkJobList = [('homo_sapiens','#2#')]
	for Species, AsId in YlkJobList:
		print('Perform Ylk for '+Species)
		Perform_MkYlk_test(Species,AsId)

	MkJobList = [('homo_sapiens','pan_troglodytes')]
	for Species,Target in MkJobList:
		print('Perform Standard Mk for '+Species+' comparing to '+Target)
		Perform_MkYlk_test(Species,Target)

