import os, glob, math, sys
import pandas as pd

def Parse_Mutation_in_VcfFile(file_name):

	flag = 1

	GenoDic = {}

	Header = ["Ref","Alt","ID","Type","Sel_coef","Position","Emerge_time","Allele_freq","Population"]

	Mut_dic = {}
	Mut_dic.setdefault('Ref',[])
	Mut_dic.setdefault('Alt',[])
	Mut_dic.setdefault('ID',[])
	Mut_dic.setdefault('Type',[])
	Mut_dic.setdefault('Sel_coef',[])
	Mut_dic.setdefault('Position',[])
	Mut_dic.setdefault('Emerge_time',[])
	Mut_dic.setdefault('Allele_freq',[])
	Mut_dic.setdefault('Population',[])

	fin = open(file_name,'r')
	for line in fin:
		if '#CHROM' in line:
			SampleList = line.strip().split('\t')[9:]
			SampleSize = len(SampleList)*2

			for IndId in SampleList:
				GenoDic.setdefault(IndId+'_1','')
				GenoDic.setdefault(IndId+'_2','')

		elif '#' in line[0]:pass
		else:

			Chrom,Pos,Id,Ref,Alt,Qual,Filt,Inf,Form = line.strip().split('\t')[:9]

			for i in range(len(Alt.split(','))):
				
				MutIds = Inf.split(';')[0].split('=')[1].split(',')
				MutType = 'm'+Inf.split(';')[5].split('=')[1].split(',')[i]
				SelCoef = float(Inf.split(';')[1].split('=')[1].split(',')[i])
				EmergeTime = int(Inf.split(';')[4].split('=')[1].split(',')[i])
				AlleleFreq = float(Inf.split(';')[6].split('=')[1].split(',')[i])/SampleSize
				Population = 'p'+Inf.split(';')[3].split('=')[1].split(',')[i]

				Mut_dic['Ref'].append(Ref)
				Mut_dic['Alt'].append(Alt.split(',')[i])
				Mut_dic['ID'].append(MutIds[i])
				Mut_dic['Type'].append(MutType)
				Mut_dic['Sel_coef'].append(SelCoef)
				Mut_dic['Position'].append(int(Pos))
				Mut_dic['Emerge_time'].append(EmergeTime)
				Mut_dic['Allele_freq'].append(AlleleFreq)
				Mut_dic['Population'].append(Population)

			Genotypes = line.strip().split('\t')[9:]
			for i in range(len(Genotypes)):
				MutState1 = int(Genotypes[i].split('|')[0])-1
				MutState2 = int(Genotypes[i].split('|')[1])-1
				if Genotypes[i].split('|')[0] != '0':GenoDic[SampleList[i]+'_1'] += MutIds[MutState1]+' '
				if Genotypes[i].split('|')[1] != '0':GenoDic[SampleList[i]+'_2'] += MutIds[MutState2]+' '
	fin.close()

	Mut_matrix = pd.DataFrame(Mut_dic, columns=Header)
	Mut_matrix = Mut_matrix.sort_values(by=['Position'])
	Mut_matrix.to_csv(file_name.split('.')[0]+'_dataframe.txt',sep='\t')

	return Mut_matrix, GenoDic


def CreateGffFile():

	GeneCoordinateDf = {};Header = ['GeneId','Start','End']
	for ColName in Header:GeneCoordinateDf.setdefault(ColName,[])

	fout = open('InitialState.gff','w');base = 0
	for i in range(1):

		GeneId = 'G'+str(i+1).zfill(4)

		### Intergenic ###
		NonCodingLen = 35000
		base += NonCodingLen

		### Genic ###
		ExonLen = 180;IntronLen = 1500
		line = '1\tDummy\tCDS\t'+\
			   str(base+1)+'\t'+\
			   str(base+ExonLen)+\
			   '\t.\t+\t.\t'+\
			   GeneId+'.fas'+'\n'
		fout.write(line)

		GeneCoordinateDf[Header[0]].append(GeneId)
		GeneCoordinateDf[Header[1]].append(base)
		GeneCoordinateDf[Header[2]].append(base+ExonLen)

		base += ExonLen

		for j in range(9):

			base += IntronLen
			line = '1\tDummy\tCDS\t'+\
				   str(base+1)+'\t'+\
				   str(base+ExonLen)+\
				   '\t.\t+\t.\t'+\
				   GeneId+'.fas'+'\n'
			fout.write(line)

			GeneCoordinateDf[Header[0]].append(GeneId)
			GeneCoordinateDf[Header[1]].append(base)
			GeneCoordinateDf[Header[2]].append(base+ExonLen)

			base += ExonLen

		### Intergenic ###
		base += NonCodingLen
	fout.close()

	GeneCoordinateDf = pd.DataFrame(GeneCoordinateDf,columns=Header)
	GeneCoordinateDf.to_csv('GeneCoDf.txt',sep='\t')


def Parse_initial_state(file_name):

	fin = open(file_name,'r')
	for line in fin:
		if 'Ancestral:' in line:
			nuc_list = line.replace('Ancestral:','').strip()
	fin.close()
	
	fout = open('InitialState.fa','w')
	fout.write('>1\n'+nuc_list)
	fout.close()


def Vcf2Fasta(FileName):
	
	OutDir = FileName.split('.')[0]

	if not OutDir in os.listdir('./'):

		RefFasta = 'InitialState.fa';RefGff = 'InitialState.gff'
		Vcf2Fasta = Vcf2FastaDir+'vcf2fasta.pl'
		Code = 'perl '+Vcf2Fasta+ \
			   ' -f '+RefFasta+ \
			   ' -v '+FileName+ \
			   ' -g '+RefGff+ \
			   ' -e CDS --phased'
		print(Code);os.system(Code)

		os.system('mv CDS '+OutDir)
	
	
if __name__ == '__main__':

	PopSize = sys.argv[1]
	Vcf2FastaDir = '../../../Program/vcf2fasta/'

	p1 = 'Generation_678000_ind_'+PopSize+'_p1_polymorphism.vcf'
	p3 = 'Generation_687000_ind_'+PopSize+'_p3_polymorphism.vcf'
	p4 = 'Generation_760000_ind_'+PopSize+'_p4_polymorphism.vcf'
	p5 = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism.vcf'
	n1 = 'Generation_0_ind_'+PopSize+'_p1_polymorphism.vcf'
	n2 = 'Generation_148000_ind_'+PopSize+'_p1_polymorphism.vcf'
	n3 = 'Generation_491000_ind_'+PopSize+'_p4_polymorphism.vcf'

	Parse_initial_state('Simulation.log')
	CreateGffFile()
	
	Flist = [p1,p3,p4,p5, n1,n2,n3]

	for FileName in Flist:
		Parse_Mutation_in_VcfFile(FileName)
		Vcf2Fasta(FileName)

