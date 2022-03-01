import sys, os, random

def ExecuteSfsCode():

	ExonNum = int(18000/ExonLen)
	PopSize = 2000
	SampSize = 100
	IterNum = 1
	Ploidy = 2
	Mu = 0.0005
	Rho = 0.0004
	
	# Fraction of variants (beneficial, deleterious)
	SelCoef_Lineage1 = (0.20,0.0)
	SelCoef_Lineage2 = (0.20,0.0)

	Args = (IterNum,PopSize,Ploidy,SampSize,Mu,Rho)
	
	PosStart = '-TW 7.4 P 3'
	PosStart2 = '-TW 7.4 P 4'
	ThreeMass1 = '5 {0} {1}'.format(SelCoef_Lineage1[0],SelCoef_Lineage1[1])
	ThreeMass2 = '5 {0} {1}'.format(SelCoef_Lineage2[0],SelCoef_Lineage2[1])

	PosSel = [PosStart+' L '+str(i)+' 1 '+ThreeMass1 for i in range(ExonNum)]
	PosSel2 = [PosStart2+' L '+str(i)+' 1 '+ThreeMass2 for i in range(ExonNum)]

	SystemStr = '../../../Program/sfscode/bin/sfs_code 7 %d '+ \
				'--length '+str(ExonNum)+' '+str(ExonLen)+' '+ \
				'--annotate F Locus.txt '+ \
				'--sex 0 '+ \
				'--popSize %d '+ \
				'--ploidy %d '+ \
				'--sampSize %d '+ \
				'--substMod 0 '+ \
				'--theta %.4f '+ \
				'--rho %.4f '+ \
				'-TS 0.0 0 1 '+ \
				'-TS 0.0 0 2 '+ \
				'-TE 0.0 0 '+ \
				'-TS 7.4 1 3 '+ \
				'-TS 7.4 1 4 '+ \
				'-TE 7.4 1 '+ \
				' '.join(PosSel)+' '+ \
				' '.join(PosSel2)+' '+ \
				'-TS 24.55 4 5 '+ \
				'-TS 24.55 4 6 '+ \
				'-TE 24.55 4 '+ \
				'-TE 33.25 5 '+ \
				'-TE 33.9 3 '+ \
				'-TE 34.35 2 '+ \
				'-TE 38 6 '+ \
				'--VCF '+ \
				'-o SfsSimulation.vcf '+ \
				'-e ErrReport.txt'

	os.system(SystemStr % Args)


def SubsetReformatVcf(VcfFile):

	fin = open(VcfFile,'r');Header = [];InitSeq = '';VcfDic = {}
	for line in fin:

		if line.startswith('##locus_'):
			InitSeq += line.strip().split('=')[1]
		elif line.startswith('##'):
			Header.append(line)
		elif line.startswith('#CHROM'):
			LastHeader = '\t'.join(line.strip().split('\t')[:9])
			for i in range(7):
				Samples = '\t'.join(line.strip().split('\t')[9+(i*N):9+(i+1)*N])
				VcfDic.setdefault(FileNameList[i],[])
				VcfDic[FileNameList[i]].append(LastHeader+'\t'+Samples+'\n')
		else:
			if line.strip() == '':continue
			Pos = int(line.split('\t')[0])*ExonLen+int(line.split('\t')[1])+1
			VarInfo1 = '1\t'+str(Pos)+'\t'+'\t'.join(line.strip().split('\t')[2:7])
			Anns = line.strip().split('\t')[7]
			Format = line.strip().split('\t')[8]

			for i in range(7):

				PopString = '\t'.join(line.strip().split('\t')[9+(i*N):9+(i+1)*N])
				NS = 'NS='+str(N)
				AF = 'AF='+str(float(PopString.count('1'))/(2*N))
				NewAnns = NS+';'+AF+';'+';'.join(Anns.split(';')[2:])
				NewLine = VarInfo1+'\t'+NewAnns+'\t'+Format+'\t'+PopString+'\n'

				if '1' in PopString:
					VcfDic[FileNameList[i]].append(NewLine)

	fin.close()

	## Create VCF files
	for FileName in FileNameList:
		fout = open(FileName,'w')
		for line in Header:
			fout.write(line)
		for VcfString in VcfDic[FileName]:
			fout.write(VcfString)
		fout.close()

	## Create Initial state sequence
	fout = open('InitialState.fa','w')
	fout.write('>1\n'+InitSeq);fout.close()

	os.system('rm SfsSimulation.vcf')


def RandomSampleVcf(VcfFile,SampleSize):

	SubVcfFile = VcfFile.replace('_100_','_'+str(SampleSize)+'_')

	fout = open(SubVcfFile,'w')
	fin = open(VcfFile,'r');Header = []
	for line in fin:

		if line.startswith('##'):
			fout.write(line)
		elif line.startswith('#CHROM'):
			TotalSamples = line.strip().split('\t')[9:]
			RanInd = random.sample(range(len(TotalSamples)),SampleSize)
			LastHeader = '\t'.join(line.strip().split('\t')[:9])
			SubSamples = [TotalSamples[i] for i in RanInd]
			fout.write(LastHeader+'\t'+'\t'.join(SubSamples)+'\n')
		else:
			VarInfo1 = '\t'.join(line.strip().split('\t')[:7])
			Anns = line.strip().split('\t')[7]
			Format = line.strip().split('\t')[8]
			Samples = line.strip().split('\t')[9:]
			SubSamples = '\t'.join([Samples[i] for i in RanInd])

			NS = 'NS='+str(SampleSize)
			AF = 'AF='+str(float(SubSamples.count('1'))/(2*SampleSize))
			NewAnns = NS+';'+AF+';'+';'.join(Anns.split(';')[2:])

			if '1' in SubSamples:
				fout.write(VarInfo1+'\t'+NewAnns+'\t'+Format+'\t'+SubSamples+'\n')

	fin.close();fout.close()


def RandomSample():

	SampleSizes = [5, 20]
	for VcfFile in FileNameList:
		for SampleSize in SampleSizes:
			RandomSampleVcf(VcfFile,SampleSize)

	
FileNameList = ['Generation_0_ind_100_p0_polymorphism.vcf',
				'Generation_148000_ind_100_p1_polymorphism.vcf',
				'Generation_687000_ind_100_p2_polymorphism.vcf',
				'Generation_678000_ind_100_p3_polymorphism.vcf',
				'Generation_491000_ind_100_p4_polymorphism.vcf',
				'Generation_665000_ind_100_p5_polymorphism.vcf',
				'Generation_760000_ind_100_p6_polymorphism.vcf']

N = 100
ExonLen = 360

ExecuteSfsCode()
SubsetReformatVcf('SfsSimulation.vcf')
RandomSample()
