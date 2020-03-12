import glob, os, multiprocessing
import pandas as pd
from scipy import stats

def ExecuteFastMlAsr(IterId,PopSize):

	## Reconstruct ancestral sequence using FastML
	## Input: Forward simulation fasta files
	## Output: '2.2.FastML_PopSize_*/' directories with ancestral sequences

	os.system('rm -r 2.2.FastML_PopSize_'+PopSize)
	os.system('mkdir 2.2.FastML_PopSize_'+PopSize)

	TreeFile = 'Simulation_tree.tre'

	RepSeqDic = {}

	P1dir = 'Generation_678000_ind_'+PopSize+'_p1_polymorphism/'
	P3dir = 'Generation_687000_ind_'+PopSize+'_p3_polymorphism/'
	P4dir = 'Generation_760000_ind_'+PopSize+'_p4_polymorphism/'
	P5dir = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism/'

	RepSeqDic['p1'] = Parse_fasta('i0_a',P1dir+'G0001.fas')
	RepSeqDic['p3'] = Parse_fasta('i0_a',P3dir+'G0001.fas')
	RepSeqDic['p4'] = Parse_fasta('i0_a',P4dir+'G0001.fas')
	RepSeqDic['p5'] = Parse_fasta('i0_a',P5dir+'G0001.fas')

	FastML_MSA_file = IterId+'_FastML_MSA.fas'
	fout = open(FastML_MSA_file,'w')
	for s_id in sorted(RepSeqDic.keys()):fout.write('>'+s_id+'\n'+RepSeqDic[s_id]+'\n')
	fout.close()

	if not FastML_MSA_file in os.listdir('./2.2.FastML_PopSize_'+PopSize+'/'):
		os.system('cp '+TreeFile+' 2.2.FastML_PopSize_'+PopSize)
		os.system('mv '+FastML_MSA_file+' 2.2.FastML_PopSize_'+PopSize)

	os.chdir('2.2.FastML_PopSize_'+PopSize)
	Commandline = 'perl ../../../'+ProgDir+'FastML.v3.11/www/fastml/FastML_Wrapper.pl'+\
				  ' --MSA_File '+FastML_MSA_file+\
				  ' --Tree Simulation_tree.tre'+\
				  ' --seqType nuc'+\
				  ' --SubMatrix JC_Nuc'+\
				  ' --indelReconstruction BOTH'+\
				  ' --outDir ./'

	print(Commandline);os.system(Commandline)

	## reduce storage
	os.system('rm prob.marginal.*')
	os.system('rm LogLikelihood_prob.margianl.csv')
	##

	os.chdir('../')


def ExecutePrankAsr(IterId,PopSize):

	## Reconstruct ancestral sequence using Prank
	## Input: Forward simulation fasta files
	## Output: '2.1.Prank_PopSize_*/' directories with ancestral sequences

	os.system('mkdir 2.1.Prank_PopSize_'+PopSize)
	
	TreeFile = 'Simulation_tree.tre'
	
	RepSeqDic = {}

	P1dir = 'Generation_678000_ind_'+PopSize+'_p1_polymorphism/'
	P3dir = 'Generation_687000_ind_'+PopSize+'_p3_polymorphism/'
	P4dir = 'Generation_760000_ind_'+PopSize+'_p4_polymorphism/'
	P5dir = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism/'

	RepSeqDic['p1'] = Parse_fasta('i0_a',P1dir+'G0001.fas')
	RepSeqDic['p3'] = Parse_fasta('i0_a',P3dir+'G0001.fas')
	RepSeqDic['p4'] = Parse_fasta('i0_a',P4dir+'G0001.fas')
	RepSeqDic['p5'] = Parse_fasta('i0_a',P5dir+'G0001.fas')
	
	PrankInput = IterId+'_Prank.fas'
	fout = open(PrankInput,'w')
	for s_id in sorted(RepSeqDic.keys()):fout.write('>'+s_id+'\n'+RepSeqDic[s_id]+'\n')
	fout.close()
	
	Commandline = '../../'+ProgDir+'prank-msa/src/prank'+\
				  ' -d='+PrankInput+\
				  ' -o=2.1.Prank_PopSize_'+PopSize+'/'+IterId+\
				  ' -f=fasta'+\
				  ' -t='+TreeFile+\
				  ' -showanc'+\
				  ' -showevents'+\
				  ' -quiet'+\
				  ' -once'+\
				  ' -kappa=1'+\
				  ' -rho=1'+\
				  ' -F'+\
				  ' -keep'

	print(Commandline);os.system(Commandline)

	os.system('rm '+PrankInput)


def ReportVariantLoci(GeneId,Sseq,Tseq,SmutDf,TmutDf):

	## Report divergence of sequence between species and polymorphism within population
	## Input: Subject and target CDS sequences and polymorphisms in dataframe
	## Output: Dataframe of variants

	MafFilter = 0.0
	b_dN = 0;b_pN = 0;dN = 0;dS = 0;pN = 0;pS = 0
	VarLociDf = {};Header = ['IterId','Pos','Sstate','Ref','Alt','SelCoef','GO','PO','VarType']
	for ColName in Header:VarLociDf.setdefault(ColName,[])

	## Filter for low minor allele frequency
	SpolDf = SmutDf[(SmutDf['Allele_freq'] < 1.0 - MafFilter) & \
					(SmutDf['Allele_freq'] > 0.0 + MafFilter)].reset_index(drop=True)

	## Count divergence and polymorphism
	for i in range(len(Sseq)):

		base = i+1
		## No Divergence ##
		if Sseq[i] == Tseq[i]:

			## Polymorphism ##
			while len(SpolDf.index) != 0 and base == SpolDf.iloc[0]['MsaPos']:

				SelCoef = float(SpolDf.iloc[0]['Sel_coef'])
				Ref = SpolDf.iloc[0]['Ref'];Alt = SpolDf.iloc[0]['Alt']

				VarLociDf[Header[0]].append(GeneId)
				VarLociDf[Header[1]].append(base)
				VarLociDf[Header[2]].append(Sseq[i])
				VarLociDf[Header[3]].append(Ref)
				VarLociDf[Header[4]].append(Alt)
				VarLociDf[Header[5]].append(SelCoef)
				VarLociDf[Header[6]].append(SpolDf.iloc[0]['Emerge_time'])
				VarLociDf[Header[7]].append(SpolDf.iloc[0]['Population'])

				if base % 4 == 0:
					pS += 1
					VarLociDf[Header[8]].append('pS')
				else:
					pN += 1
					VarLociDf[Header[8]].append('pN')
					if SelCoef > 0:b_pN += 1

				SpolDf = SpolDf.drop([0]).reset_index(drop=True)

		## Divergence ##
		elif Sseq[i] != Tseq[i]:
			Smut = SmutDf[(SmutDf['MsaPos'] == base) & (SmutDf['Alt'] == Sseq[i])]

			if len(SpolDf.index) != 0:

				## No Polymorphism ##
				if base != SpolDf.iloc[0]['MsaPos']:

					VarLociDf[Header[0]].append(GeneId)
					VarLociDf[Header[1]].append(base)
					VarLociDf[Header[2]].append(Sseq[i])
					VarLociDf[Header[3]].append(Sseq[i])
					VarLociDf[Header[4]].append(Tseq[i])

					SelCoef = 0;GO = 0;PO='p1'
					if len(Smut.index) != 0:
						SelCoef = float(Smut['Sel_coef'].values[0])
						GO = int(Smut['Emerge_time'].values[0])
						PO = str(Smut['Population'].values[0])

					VarLociDf[Header[5]].append(SelCoef)
					VarLociDf[Header[6]].append(GO)
					VarLociDf[Header[7]].append(PO)
					if base % 4 == 0:
						dS += 1
						VarLociDf[Header[8]].append('dS')
					else:
						dN += 1
						VarLociDf[Header[8]].append('dN')
						if SelCoef > 0:b_dN += 1

				## Polymorphism ##
				while len(SpolDf.index) != 0 and base == SpolDf.iloc[0]['MsaPos']:

					SelCoef = float(SpolDf.iloc[0]['Sel_coef'])
					Ref = SpolDf.iloc[0]['Ref'];Alt = SpolDf.iloc[0]['Alt']

					VarLociDf[Header[0]].append(GeneId)
					VarLociDf[Header[1]].append(base)
					VarLociDf[Header[2]].append(Sseq[i])
					VarLociDf[Header[3]].append(Ref)
					VarLociDf[Header[4]].append(Alt)
					VarLociDf[Header[5]].append(SelCoef)
					VarLociDf[Header[6]].append(SpolDf.iloc[0]['Emerge_time'])
					VarLociDf[Header[7]].append(SpolDf.iloc[0]['Population'])

					if base % 4 == 0:
						pS += 1
						VarLociDf[Header[8]].append('pS')
					else:
						pN += 1
						VarLociDf[Header[8]].append('pN')
						if SelCoef > 0:b_pN += 1

					if not Tseq[i] in [Ref,Alt]:

						VarLociDf[Header[0]].append(GeneId)
						VarLociDf[Header[1]].append(base)
						VarLociDf[Header[2]].append(Sseq[i])
						VarLociDf[Header[3]].append(Sseq[i])
						VarLociDf[Header[4]].append(Tseq[i])

						SelCoef = 0;GO = 0;PO='p1'
						if len(Smut.index) != 0:
							SelCoef = float(Smut['Sel_coef'].values[0])
							GO = int(Smut['Emerge_time'].values[0])
							PO = str(Smut['Population'].values[0])

						VarLociDf[Header[5]].append(SelCoef)
						VarLociDf[Header[6]].append(GO)
						VarLociDf[Header[7]].append(PO)

						if base % 4 == 0:
							dS += 1
							VarLociDf[Header[8]].append('dS')
						else:
							dN += 1
							VarLociDf[Header[8]].append('dN')
							if SelCoef > 0:b_dN += 1

					SpolDf = SpolDf.drop([0]).reset_index(drop=True)

			## No polymorphism ##
			elif len(SpolDf.index) == 0:

				VarLociDf[Header[0]].append(GeneId)
				VarLociDf[Header[1]].append(base)
				VarLociDf[Header[2]].append(Sseq[i])
				VarLociDf[Header[3]].append(Sseq[i])
				VarLociDf[Header[4]].append(Tseq[i])

				SelCoef = 0;GO = 0;PO='p1'
				if len(Smut.index) != 0:
					SelCoef = float(Smut['Sel_coef'].values[0])
					GO = int(Smut['Emerge_time'].values[0])
					PO = str(Smut['Population'].values[0])

				VarLociDf[Header[5]].append(SelCoef)
				VarLociDf[Header[6]].append(GO)
				VarLociDf[Header[7]].append(PO)
				if base % 4 == 0:
					dS += 1
					VarLociDf[Header[8]].append('dS')
				else:
					dN += 1
					VarLociDf[Header[8]].append('dN')
					if SelCoef > 0:b_dN += 1

	VarLociDf = pd.DataFrame(VarLociDf,columns=Header)

	return VarLociDf, dN,dS,pN,pS,b_dN,b_pN


def Parse_reconstructed_seq(fname):

	## Parse Fasta file containing ancestral sequence

	seq_dic = {}
	fin = open(fname,'r')
	for line in fin:
		if '>' in line[0]:
			fsid = line.strip().split(' ')[0].replace('>','')
			seq_dic.setdefault(fsid,'')
		else:seq_dic[fsid] += line.strip()
	fin.close()

	return seq_dic


def Parse_fasta(s_id, fasta_file):

	## Parse specific sequence ID from a Fasta file

	sequence = ''

	flag = 0
	fin = open(fasta_file,'r')
	for line in fin:
		if '>' in line[0]:
			c_s_id = line.strip().replace('>','')
			if c_s_id == s_id:flag = 1
			else:flag = 0
		else:
			if flag == 1: sequence += line.strip()
	fin.close()

	return sequence


def ModifyMutDf(PopSize):

	## Add gene info into polymorphism dataframe

	flist = sorted(glob.glob('*_ind_'+PopSize+'_*_polymorphism_dataframe.txt'))
	for f in flist:

		CumulativeLenDic = {};VarMsaPosList = [];VarGeneList = []
		GeneCoDf = pd.read_csv('GeneCoDf.txt',sep='\t',index_col=0)
		CurGene = GeneCoDf.iloc[0]['GeneId']
		CurStart = GeneCoDf.iloc[0]['Start']
		CurEnd = GeneCoDf.iloc[0]['End']

		CumulativeLenDic.setdefault(CurGene,0)

		print('Modyfying %s' % f)

		MutDf = pd.read_csv(f,sep='\t',index_col=0)
		MutDf = MutDf.sort_values(by=['Position']).reset_index(drop=True)

		for i in range(len(MutDf.index)):

			Pos = MutDf.iloc[i]['Position'] - 1

			while Pos > CurEnd and len(GeneCoDf.index) != 0:

				CumulativeLenDic[CurGene] += CurEnd - CurStart
				GeneCoDf = GeneCoDf.drop([0]).reset_index(drop=True)
				if len(GeneCoDf.index) != 0:

					CurGene = GeneCoDf.iloc[0]['GeneId']
					CurStart = GeneCoDf.iloc[0]['Start']
					CurEnd = GeneCoDf.iloc[0]['End']

					CumulativeLenDic.setdefault(CurGene,0)

			if Pos >= CurStart and Pos < CurEnd:

				NewPos = CumulativeLenDic[CurGene] + Pos - CurStart + 1
				VarMsaPosList.append(NewPos)
				VarGeneList.append(CurGene)

		MutDf['MsaPos'] = VarMsaPosList
		MutDf['IterId'] = VarGeneList

		MutDf.to_csv(f.replace('_dataframe.txt','_ModDataframe.txt'),sep='\t')


def PrankMkYlk(IterId,Subject,Target, PopSize):

	## Count divergences (dN, dS) and polymorphisms (pN, pS) and calculate MK score
	## Input: Variant calling files and Fasta file containing aligned sequences
	## Output: Variant report dataframe and Mk/Ylk dataframe

	IdDic = {'p1':'p1', 'p3':'p3', 'p4':'p4', 'p5':'p5', '#3#':'1', '#2#':'2', '#1#':'3'}

	p1 = 'Generation_678000_ind_'+PopSize+'_p1_polymorphism.vcf'
	p3 = 'Generation_687000_ind_'+PopSize+'_p3_polymorphism.vcf'
	p4 = 'Generation_760000_ind_'+PopSize+'_p4_polymorphism.vcf'
	p5 = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism.vcf'

	n1 = 'Generation_0_ind_'+PopSize+'_p1_polymorphism.vcf'
	n2 = 'Generation_148000_ind_'+PopSize+'_p1_polymorphism.vcf'
	n3 = 'Generation_491000_ind_'+PopSize+'_p4_polymorphism.vcf'

	PopDic = {'p1':p1, 'p3':p3, 'p4':p4, 'p5':p5, '#3#':n1, '#2#':n2, '#1#':n3}

	SubjectVcf = PopDic[Subject];TargetVcf = PopDic[Target]
	SmutDf = pd.read_csv(SubjectVcf.replace('.vcf','_ModDataframe.txt'),sep='\t',index_col=0)
	TmutDf = pd.read_csv(TargetVcf.replace('.vcf','_ModDataframe.txt'),sep='\t',index_col=0)

	FileName = '2.1.Prank_PopSize_'+PopSize+'/'+IterId+'.anc.fas'
	SeqDic = Parse_reconstructed_seq(FileName)
	Sseq = SeqDic[Subject]
	Tseq = SeqDic[Target]

	MktDf = {};MktHeader = ['IterId','dN','dS','pN','pS','bdN','bpN','MkScore','FSPVal']
	for ColName in MktHeader:MktDf.setdefault(ColName,[])

	GeneVarDf,dN,dS,pN,pS,bdN,bpN = ReportVariantLoci(IterId,Sseq,Tseq,SmutDf,TmutDf)

	#if any(x < 5 for x in [dN,dS,pN,pS]):
	if 0 in [dN,dS,pN,pS]:
		PseudoCount = 1.0
	else:
		PseudoCount = 0.0
	MkScore = ((dN+PseudoCount)/(dS+PseudoCount))/((pN+PseudoCount)/(pS+PseudoCount))
	OddRatio, FSPVal = stats.fisher_exact([[dN, dS], [pN, pS]])

	MktDf[MktHeader[0]].append(IterId)
	MktDf[MktHeader[1]].append(dN)
	MktDf[MktHeader[2]].append(dS)
	MktDf[MktHeader[3]].append(pN)
	MktDf[MktHeader[4]].append(pS)
	MktDf[MktHeader[5]].append(bdN)
	MktDf[MktHeader[6]].append(bpN)
	MktDf[MktHeader[7]].append(MkScore)
	MktDf[MktHeader[8]].append(FSPVal)

	Subject = IdDic[Subject];Target = IdDic[Target]
	MktDf = pd.DataFrame(MktDf,columns=MktHeader)
	MktDf.to_csv(Subject+'_vs_'+Target+'_Popsize'+PopSize+'_Prank-MkYlk.txt',sep='\t')
	GeneVarDf.to_csv(Subject+'_vs_'+Target+'_Popsize'+PopSize+'_Prank-VarRept.txt',sep='\t')


def FastMlMkYlk(IterId,Subject,Target, PopSize):

	## Count divergences (dN, dS) and polymorphisms (pN, pS) and calculate MK score
	## Input: Variant calling files and Fasta file containing aligned sequences
	## Output: Variant report dataframe and Mk/Ylk dataframe

	p1 = 'Generation_678000_ind_'+PopSize+'_p1_polymorphism.vcf'
	p3 = 'Generation_687000_ind_'+PopSize+'_p3_polymorphism.vcf'
	p4 = 'Generation_760000_ind_'+PopSize+'_p4_polymorphism.vcf'
	p5 = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism.vcf'

	n1 = 'Generation_0_ind_'+PopSize+'_p1_polymorphism.vcf'
	n2 = 'Generation_148000_ind_'+PopSize+'_p1_polymorphism.vcf'
	n3 = 'Generation_491000_ind_'+PopSize+'_p4_polymorphism.vcf'

	PopDic = {'p1':p1, 'p3':p3, 'p4':p4, 'p5':p5, '1':n1, '2':n2, '3':n3}

	SubjectVcf = PopDic[Subject];TargetVcf = PopDic[Target]
	SmutDf = pd.read_csv(SubjectVcf.replace('.vcf','_ModDataframe.txt'),sep='\t',index_col=0)
	TmutDf = pd.read_csv(TargetVcf.replace('.vcf','_ModDataframe.txt'),sep='\t',index_col=0)

	FileName = '2.2.FastML_PopSize_'+PopSize+'/FilesForJalView/seq.joint_JalView.FASTA.aln'
	SeqDic = Parse_reconstructed_seq(FileName)
	Sseq = SeqDic[Subject]
	Tseq = SeqDic[Target]

	MktDf = {};MktHeader = ['IterId','dN','dS','pN','pS','bdN','bpN','MkScore','FSPVal']
	for ColName in MktHeader:MktDf.setdefault(ColName,[])

	GeneVarDf,dN,dS,pN,pS,bdN,bpN = ReportVariantLoci(IterId,Sseq,Tseq,SmutDf,TmutDf)

	#if any(x < 5 for x in [dN,dS,pN,pS]):
	if 0 in [dN,dS,pN,pS]:
		PseudoCount = 1.0
	else:
		PseudoCount = 0.0
	MkScore = ((dN+PseudoCount)/(dS+PseudoCount))/((pN+PseudoCount)/(pS+PseudoCount))
	OddRatio, FSPVal = stats.fisher_exact([[dN, dS], [pN, pS]])

	MktDf[MktHeader[0]].append(IterId)
	MktDf[MktHeader[1]].append(dN)
	MktDf[MktHeader[2]].append(dS)
	MktDf[MktHeader[3]].append(pN)
	MktDf[MktHeader[4]].append(pS)
	MktDf[MktHeader[5]].append(bdN)
	MktDf[MktHeader[6]].append(bpN)
	MktDf[MktHeader[7]].append(MkScore)
	MktDf[MktHeader[8]].append(FSPVal)

	MktDf = pd.DataFrame(MktDf,columns=MktHeader)
	MktDf.to_csv(Subject+'_vs_'+Target+'_Popsize'+PopSize+'_FastMl-MkYlk.txt',sep='\t')
	GeneVarDf.to_csv(Subject+'_vs_'+Target+'_Popsize'+PopSize+'_Fastml-VarRept.txt',sep='\t')


def ExecuteMultipleSimulation(Directory,SimNum,ThrNum):

	## Perform Simulation in parallel
	## Simulation options are prepared as in 'TemplateData/' Directory
	## Input: Simulation Directory containing 'TemplateData/' Directory
	## Output: Simulation data representing population of individuals their genes

	os.chdir(Directory);IterDirList = []

	for i in range(SimNum):

		IterId = 'Iter'+str(i+1).zfill(4)
		IterDirList.append(Directory+'/'+IterId)

		if not IterId in os.listdir('./'):os.system('cp -R TemplateData '+IterId)

	p = multiprocessing.Pool(ThrNum)
	p.map(ExecuteSimulation,IterDirList)

	os.chdir('../')


def ExecuteSimulation(Directory):

	## Perform Simulation using the option described in 'Chrom_sim_shorter_ver3.eidos' file
	## Input: Copies of 'TemplateData/' directory containing simulation eidos files
	## Output: Simulation data and Mk/Ylk analysis results

	os.chdir(Directory.split('/')[-1])

	RunId = Directory.split('/')[-2]+' '+Directory.split('/')[-1]
	systemstr = '../../'+ProgDir+'SLiM/build/slim ./Chrom_sim_shorter_ver3.eidos > Simulation.log'
	#print(RunId+' Slim Simulation');os.system(systemstr)
	
	PopSizeList = ['5','20','100','500']
	for PopSize in PopSizeList:

		systemstr = 'python 001.Output_fasta_file_ver4.py '+PopSize
		print(RunId+' Create Fasta File');os.system(systemstr)
	
		ModifyMutDf(PopSize)
		ExecutePrankAsr(Directory.split('/')[-1],PopSize)
		ExecuteFastMlAsr(Directory.split('/')[-1],PopSize) #

		## Execute YLK - Prank
		PrankMkYlk(Directory.split('/')[-1],'p1','#2#',PopSize)
		## Execute YLK - FastMl
		FastMlMkYlk(Directory.split('/')[-1],'p1','2',PopSize) #
		## Execute standard MK
		FastMlMkYlk(Directory.split('/')[-1],'p1','p5',PopSize) #

	os.chdir('../')


def MakePrediction(MkYlkDf,SelThr):

	Divergence = 34800
	SelCoefList = [];PosSelList = [];TP = 0;TN = 0;FP = 0;FN = 0

	p1 = 'Generation_678000_ind_'+PopSize+'_p1_polymorphism_ModDataframe.txt'
	p3 = 'Generation_687000_ind_'+PopSize+'_p3_polymorphism_ModDataframe.txt'
	p4 = 'Generation_760000_ind_'+PopSize+'_p4_polymorphism_ModDataframe.txt'
	p5 = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism_ModDataframe.txt'

	n1 = 'Generation_0_ind_'+PopSize+'_p1_polymorphism_ModDataframe.txt'
	n2 = 'Generation_148000_ind_'+PopSize+'_p1_polymorphism_ModDataframe.txt'
	n3 = 'Generation_491000_ind_'+PopSize+'_p4_polymorphism_ModDataframe.txt'

	PopDic = {'p1':p1, 'p3':p3, 'p4':p4, 'p5':p5, '1':n1, '2':n2, '3':n3}

	for i in MkYlkDf.index:

		IterId = MkYlkDf.loc[i,'IterId']
		MktScore = MkYlkDf.loc[i,'MkScore']
		VarDfFile = PopDic['p1']
		VarDf = pd.read_csv(IterId+'/'+VarDfFile,sep='\t',index_col=0)
		RedVarDf = VarDf[(VarDf['Emerge_time'] >= Divergence) & (VarDf['Allele_freq'] ==1.0)]

		if sum(RedVarDf['Sel_coef']) > SelThr:
			PosSel = 1
			if MktScore > 1:
				TP += 1
			else:
				FN += 1
		else:
			PosSel = 0
			if MktScore > 1:
				FP += 1
			else:
				TN += 1

		SelCoefList.append(sum(RedVarDf['Sel_coef']))
		PosSelList.append(PosSel)

	return TP,FN,FP,TN,SelCoefList,PosSelList


def SummarizeAllSimulation(PopSize,Directory):

	## Summarize the overall results of simulation analysis
	## Input: Simulation data and Mk/Ylk analysis results
	## Output: Performance of Mk/Ylk analyses

	os.chdir(Directory)

	PredictionDf = {}
	Header = ['Method','Selection','PopSize','TP','TN','FP','FN','Sen','Spe','Acc']
	for ColName in Header:PredictionDf.setdefault(ColName,[])

	IterList = sorted(glob.glob('Iter*'));SelThrList = [0]

	for SelThr in SelThrList:

		## Novel MK ##

		FastMlYlkDf = pd.DataFrame();PrankYlkDf = pd.DataFrame()
		for Element in IterList:

			IterId = Element.split('/')[-1]
			MktFile = IterId+'/p1_vs_2_Popsize'+PopSize+'_FastMl-MkYlk.txt'
			SubDf = pd.read_csv(MktFile,sep='\t',index_col=0)
			FastMlYlkDf = FastMlYlkDf.append(SubDf)

			MktFile = IterId+'/p1_vs_2_Popsize'+PopSize+'_Prank-MkYlk.txt'
			SubDf = pd.read_csv(MktFile,sep='\t',index_col=0)
			PrankYlkDf = PrankYlkDf.append(SubDf)

		## FastMl Ylk	
		FastMlYlkDf = FastMlYlkDf.reset_index(drop=True)
		TP,FN,FP,TN,SelCoefList,PosSelList = MakePrediction(FastMlYlkDf,SelThr)

		Sen = 'NA';Spe = 'NA';Acc = 'NA'
		if not TP+FN == 0:
			Sen = float(TP)/float(TP+FN)
		if not TN+FP == 0:
			Spe = float(TN)/float(TN+FP)
		Acc = float(TP+TN)/float(TP+TN+FP+FN)

		FastMlYlkDf['SelCoef'] = SelCoefList
		FastMlYlkDf['PosSel'] = PosSelList

		OutFile = 'p1_vs_2_SelThr'+str(SelThr)+'_PopSize'+str(PopSize)+'_FastMl-MkYlk-f.txt'
		FastMlYlkDf.to_csv(OutFile,sep='\t')

		PredictionDf[Header[0]].append('YLK-FastMl')
		PredictionDf[Header[1]].append(Directory)
		PredictionDf[Header[2]].append(PopSize)
		PredictionDf[Header[3]].append(TP)
		PredictionDf[Header[4]].append(TN)
		PredictionDf[Header[5]].append(FP)
		PredictionDf[Header[6]].append(FN)
		PredictionDf[Header[7]].append(Sen)
		PredictionDf[Header[8]].append(Spe)
		PredictionDf[Header[9]].append(Acc)

		## Prank Ylk
		PrankYlkDf = PrankYlkDf.reset_index(drop=True)
		TP,FN,FP,TN,SelCoefList,PosSelList = MakePrediction(PrankYlkDf,SelThr)

		Sen = 'NA';Spe = 'NA';Acc = 'NA'
		if not TP+FN == 0:
			Sen = float(TP)/float(TP+FN)
		if not TN+FP == 0:
			Spe = float(TN)/float(TN+FP)
		Acc = float(TP+TN)/float(TP+TN+FP+FN)

		PrankYlkDf['SelCoef'] = SelCoefList
		PrankYlkDf['PosSel'] = PosSelList

		OutFile = 'p1_vs_2_SelThr'+str(SelThr)+'_PopSize'+str(PopSize)+'_Prank-MkYlk-f.txt'
		PrankYlkDf.to_csv(OutFile,sep='\t')

		PredictionDf[Header[0]].append('YLK-Prank')
		PredictionDf[Header[1]].append(Directory)
		PredictionDf[Header[2]].append(PopSize)
		PredictionDf[Header[3]].append(TP)
		PredictionDf[Header[4]].append(TN)
		PredictionDf[Header[5]].append(FP)
		PredictionDf[Header[6]].append(FN)
		PredictionDf[Header[7]].append(Sen)
		PredictionDf[Header[8]].append(Spe)
		PredictionDf[Header[9]].append(Acc)

		## Standard MK ##

		StdMktDf = pd.DataFrame()
		for IterId in IterList:
		
			MktFile = IterId+'/p1_vs_p5_Popsize'+PopSize+'_FastMl-MkYlk.txt'
			SubDf = pd.read_csv(MktFile,sep='\t',index_col=0)
			StdMktDf = StdMktDf.append(SubDf)
	
		StdMktDf = StdMktDf.reset_index(drop=True)
		TP,FN,FP,TN,SelCoefList,PosSelList = MakePrediction(StdMktDf,SelThr)

		Sen = 'NA';Spe = 'NA';Acc = 'NA'
		if not TP+FN == 0:
			Sen = float(TP)/float(TP+FN)
		if not TN+FP == 0:
			Spe = float(TN)/float(TN+FP)
		Acc = float(TP+TN)/float(TP+TN+FP+FN)

		StdMktDf['SelCoef'] = SelCoefList
		StdMktDf['PosSel'] = PosSelList

		OutFile = 'p1_vs_p5_SelThr'+str(SelThr)+'_PopSize'+str(PopSize)+'_FastMl-MkYlk-f.txt'
		StdMktDf.to_csv(OutFile,sep='\t')

		PredictionDf[Header[0]].append('MK')
		PredictionDf[Header[1]].append(Directory)
		PredictionDf[Header[2]].append(PopSize)
		PredictionDf[Header[3]].append(TP)
		PredictionDf[Header[4]].append(TN)
		PredictionDf[Header[5]].append(FP)
		PredictionDf[Header[6]].append(FN)
		PredictionDf[Header[7]].append(Sen)
		PredictionDf[Header[8]].append(Spe)
		PredictionDf[Header[9]].append(Acc)

	PredictionDf = pd.DataFrame(PredictionDf,columns=Header)

	os.chdir('../')

	return PredictionDf


if __name__ == '__main__':

	DirectoryList = sorted(glob.glob('Scenario-*'))

	MainDir = './'
	ProgDir = '../Program/'
	SimNum = 1000
	ThrNum = 10

	BigPredictionDf = pd.DataFrame()

	for Directory in DirectoryList:

		ExecuteMultipleSimulation(Directory,SimNum,ThrNum)

		PopSizeList = ['5','20','100','500']
		for PopSize in PopSizeList:

			PredictionDf = SummarizeAllSimulation(PopSize,Directory)
			BigPredictionDf = BigPredictionDf.append(PredictionDf)

	BigPredictionDf = BigPredictionDf.reset_index(drop=True)
	BigPredictionDf.to_csv('AllSimulations_PredictionSummary.txt',sep='\t')
