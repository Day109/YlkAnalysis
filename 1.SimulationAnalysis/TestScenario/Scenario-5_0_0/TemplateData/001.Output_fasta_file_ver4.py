import os, glob, math, sys
import pandas as pd

def Parse_Mutation_in_VcfFile(file_name):

	flag = 1

	Header = ["Ref","Alt","ID","Type","Sel_coef","Position","Emerge_time","Allele_freq"]

	Mut_dic = {}
	Mut_dic.setdefault('Ref',[])
	Mut_dic.setdefault('Alt',[])
	Mut_dic.setdefault('ID',[])
	Mut_dic.setdefault('Type',[])
	Mut_dic.setdefault('Sel_coef',[])
	Mut_dic.setdefault('Position',[])
	Mut_dic.setdefault('Emerge_time',[])
	Mut_dic.setdefault('Allele_freq',[])

	fin = open(file_name,'r')
	for line in fin:
		if '#CHROM' in line:
			SampleList = line.strip().split('\t')[9:]
			SampleSize = len(SampleList)*2

		elif '#' in line[0]:pass
		else:

			Chrom,Pos,Id,Ref,Alt,Qual,Filt,Inf,Form = line.strip().split('\t')[:9]

			MutType = Inf.split(';NSYN=')[1].split(';')[0]
			SelCoef = float(Inf.split(';S=')[1])
			EmergeTime = int(Inf.split(';GA=')[1].split(';')[0])
			MutId = Pos+'_'+Inf.split(';GA=')[1].split(';')[0]
			AlleleFreq = float(Inf.split(';AF=')[1].split(';')[0])

			Mut_dic['Ref'].append(Ref)
			Mut_dic['Alt'].append(Alt)
			Mut_dic['ID'].append(MutId)
			Mut_dic['Type'].append(MutType)
			Mut_dic['Sel_coef'].append(SelCoef)
			Mut_dic['Position'].append(int(Pos))
			Mut_dic['Emerge_time'].append(EmergeTime)
			Mut_dic['Allele_freq'].append(AlleleFreq)

	fin.close()

	Mut_matrix = pd.DataFrame(Mut_dic, columns=Header)
	Mut_matrix = Mut_matrix.sort_values(by=['Position'])
	Mut_matrix.to_csv(file_name.split('.')[0]+'_dataframe.txt',sep='\t')

	return Mut_matrix


def CreateGffFile():

	GeneCoordinateDf = {};Header = ['GeneId','Start','End']
	for ColName in Header:GeneCoordinateDf.setdefault(ColName,[])

	fout = open('InitialState.gff','w');base = 0
	for i in range(10):

		GeneId = 'G'+str(1).zfill(4)

		### Intergenic ###
		NonCodingLen = 0
		base += NonCodingLen

		### Genic ###
		ExonLen = 360;IntronLen = 0
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

		for j in range(4):

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


def Vcf2Fasta(FileName):
	
	OutDir = FileName.split('.')[0]

	if not OutDir in os.listdir('./'):

		RefFasta = 'InitialState.fa';RefGff = 'InitialState.gff'
		Vcf2Fasta = Vcf2FastaDir+'vcf2fasta.pl'
		Code = 'perl '+Vcf2Fasta+ \
			   ' -f '+RefFasta+ \
			   ' -v '+FileName+ \
			   ' -g '+RefGff+ \
			   ' -e CDS --phased > Vcf2Fasta.log 2> Vcf2Fasta.log2'
		os.system(Code)

		os.system('mv CDS '+OutDir)
	
	
if __name__ == '__main__':

	PopSize = sys.argv[1]
	Vcf2FastaDir = '../../../Program/vcf2fasta/'

	p3 = 'Generation_678000_ind_'+PopSize+'_p3_polymorphism.vcf'
	p2 = 'Generation_687000_ind_'+PopSize+'_p2_polymorphism.vcf'
	p6 = 'Generation_760000_ind_'+PopSize+'_p6_polymorphism.vcf'
	p5 = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism.vcf'
	p0 = 'Generation_0_ind_'+PopSize+'_p0_polymorphism.vcf'
	p1 = 'Generation_148000_ind_'+PopSize+'_p1_polymorphism.vcf'
	p4 = 'Generation_491000_ind_'+PopSize+'_p4_polymorphism.vcf'

	CreateGffFile()
	
	Flist = [p3,p2,p6,p5, p0,p1,p4]

	for FileName in Flist:
		Parse_Mutation_in_VcfFile(FileName)
		Vcf2Fasta(FileName)

