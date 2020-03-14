import glob,os
import pandas as pd

def CountPosGenes(DirectorySet):

	P1_N = 0;P5_N = 0;P1P5_N = 0;Non_N = 0
	P1_S = 0;P5_S = 0;P1P5_S = 0;Non_S = 0
	P1 = 0;P5 = 0;P1P5 = 0;Non = 0;Cnt = 0
	for SubDir in DirectorySet:

		IterList = sorted(glob.glob(SubDir+'/Iter*'))
		for IterDir in IterList:
			
			Cnt += 1
			if Cnt % 100 == 0:print(str(Cnt)+'th gene counted')

			p1VarFile = IterDir+'/Generation_678000_ind_100_p1_polymorphism_ModDataframe.txt'
			p5VarFile = IterDir+'/Generation_665000_ind_100_p5_polymorphism_ModDataframe.txt'
			
			p1VarDf = pd.read_csv(p1VarFile,sep='\t',index_col=0)
			Redp1VarDf = p1VarDf[(p1VarDf['Emerge_time'] >= Divergence) & (p1VarDf['Allele_freq']==1.0)]
			p5VarDf = pd.read_csv(p5VarFile,sep='\t',index_col=0)
			Redp5VarDf = p5VarDf[(p5VarDf['Emerge_time'] >= Divergence) & (p5VarDf['Allele_freq']==1.0)]

			p1FitEff = sum(Redp1VarDf['Sel_coef'])
			p5FitEff = sum(Redp5VarDf['Sel_coef'])

			NmktFile = IterDir+'/p1_vs_2_Popsize100_FastMl-MkYlk.txt'
			NmktDf = pd.read_csv(NmktFile,sep='\t',index_col=0)
			Nmkt = NmktDf.iloc[0]['MkScore']

			SmktFile = IterDir+'/p1_vs_p5_Popsize100_FastMl-MkYlk.txt'
			SmktDf = pd.read_csv(SmktFile,sep='\t',index_col=0)
			Smkt = SmktDf.iloc[0]['MkScore']

			if p1FitEff > 0:
				if p5FitEff > 0:
					P1P5 += 1
					if Nmkt > 1:P1P5_N += 1
					if Smkt > 1:P1P5_S += 1
				else:
					P1 += 1
					if Nmkt > 1:P1_N += 1
					if Smkt > 1:P1_S += 1
			else:
				if p5FitEff > 0:
					P5 += 1
					if Nmkt <= 1:P5_N += 1
					if Smkt <= 1:P5_S += 1
				else:
					Non += 1
					if Nmkt <= 1:Non_N += 1
					if Smkt <= 1:Non_S += 1

	return P1,P5,P1P5,Non, P1_N,P5_N,P1P5_N,Non_N, P1_S,P5_S,P1P5_S,Non_S


if __name__ == '__main__':

	os.chdir('../')

	CountDf = {};Header = ['Rate','Selection','P1','P5','P1P5','None','P1_N','P5_N','P1P5_N','None_N','P1_S','P5_S','P1P5_S','Non_S']
	for ColName in Header:CountDf.setdefault(ColName,[])

	RateList = glob.glob('MutRate_*')
	for Rate in RateList:

		Dir = Rate;os.chdir(Dir)
		Divergence = 34800

		DirectorySetList = [['Scenario-1_None','Scenario-2_P1','Scenario-3_P5','Scenario-4_P1P5']]

		for DirectorySet in DirectorySetList:
			
			print('Parsing genes in Simulation: '+Rate)
			P1,P5,P1P5,Non, P1_N,P5_N,P1P5_N,Non_N, P1_S,P5_S,P1P5_S,Non_S = CountPosGenes(DirectorySet)

			CountDf[Header[0]].append(Rate)
			CountDf[Header[1]].append(Rate)
			CountDf[Header[2]].append(P1)
			CountDf[Header[3]].append(P5)
			CountDf[Header[4]].append(P1P5)
			CountDf[Header[5]].append(Non)

			CountDf[Header[6]].append(P1_N)
			CountDf[Header[7]].append(P5_N)
			CountDf[Header[8]].append(P1P5_N)
			CountDf[Header[9]].append(Non_N)
			CountDf[Header[10]].append(P1_S)
			CountDf[Header[11]].append(P5_S)
			CountDf[Header[12]].append(P1P5_S)
			CountDf[Header[13]].append(Non_S)

		os.chdir('../')

	CountDf = pd.DataFrame(CountDf,columns=Header)
	CountDf.to_csv('OtherScripts/PosGeneCount_CombinedModel.txt',sep='\t')

