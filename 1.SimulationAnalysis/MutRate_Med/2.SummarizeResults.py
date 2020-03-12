import os, glob
import pandas as pd

def SummarizeAllSimulation(PopSize,DirectorySet):

	## Summarize the merged simulation results containing different selection scenarios
	## Input: Simulation data under different selection scenarios
	## Output: Overall performance summarized in dataframe

	JobId = '-'.join([i.split('_')[-1] for i in DirectorySet])

	PredictionDf = {};Header = ['Comparison','TP','TN','FP','FN','Sen','Spe','Acc']
	for ColName in Header:PredictionDf.setdefault(ColName,[])

	### Novel MK ###

	print('Parsing MKDfs for Novel MK')
	NovMktDf = pd.DataFrame()
	for Directory in DirectorySet:

		MktFile = Directory+'/p1_vs_2_SelThr0_PopSize'+PopSize+'_FastMl-MkYlk-f.txt'
		SubDf = pd.read_csv(MktFile,sep='\t',index_col=0)
		NovMktDf = NovMktDf.append(SubDf)

	NovMktDf = NovMktDf.reset_index(drop=True)

	SelThrList = [0]
	for SelThr in SelThrList:

		TP = len(NovMktDf[(NovMktDf['SelCoef'] > SelThr) & (NovMktDf['MkScore'] > 1)])
		FN = len(NovMktDf[(NovMktDf['SelCoef'] > SelThr) & (NovMktDf['MkScore'] <= 1)])
		FP = len(NovMktDf[(NovMktDf['SelCoef'] <= SelThr) & (NovMktDf['MkScore'] > 1)])
		TN = len(NovMktDf[(NovMktDf['SelCoef'] <= SelThr) & (NovMktDf['MkScore'] <= 1)])

		Sen = 'NaN';Spe = 'NaN';Acc = 'NaN'
		if not TP+FN == 0:
			Sen = float(TP)/float(TP+FN)
		if not TN+FP == 0:
			Spe = float(TN)/float(TN+FP)
		Acc = float(TP+TN)/float(TP+TN+FP+FN)

		OutFile = JobId+'_p1_vs_2_PopSize'+str(PopSize)+'_FastMl-MkYlk-f.txt'
		NovMktDf.to_csv(OutFile,sep='\t')

		Comparison = JobId+'_p1_vs_2_SelThr'+str(SelThr)+'_PopSize'+str(PopSize)
		PredictionDf[Header[0]].append(Comparison)
		PredictionDf[Header[1]].append(TP)
		PredictionDf[Header[2]].append(TN)
		PredictionDf[Header[3]].append(FP)
		PredictionDf[Header[4]].append(FN)
		PredictionDf[Header[5]].append(Sen)
		PredictionDf[Header[6]].append(Spe)
		PredictionDf[Header[7]].append(Acc)

		### Standard MK ###

	print('Parsing MKDfs for Standard MK')
	StdMktDf = pd.DataFrame()
	for Directory in DirectorySet:

		MktFile = Directory+'/p1_vs_p5_SelThr0_PopSize'+PopSize+'_FastMl-MkYlk-f.txt'
		SubDf = pd.read_csv(MktFile,sep='\t',index_col=0)
		StdMktDf = StdMktDf.append(SubDf)

	StdMktDf = StdMktDf.reset_index(drop=True)

	for SelThr in SelThrList:

		TP = len(StdMktDf[(StdMktDf['SelCoef'] > SelThr) & (StdMktDf['MkScore'] > 1)])
		FN = len(StdMktDf[(StdMktDf['SelCoef'] > SelThr) & (StdMktDf['MkScore'] <= 1)])
		FP = len(StdMktDf[(StdMktDf['SelCoef'] <= SelThr) & (StdMktDf['MkScore'] > 1)])
		TN = len(StdMktDf[(StdMktDf['SelCoef'] <= SelThr) & (StdMktDf['MkScore'] <= 1)])

		Sen = 'NaN';Spe = 'NaN';Acc = 'NaN'
		if not TP+FN == 0:
			Sen = float(TP)/float(TP+FN)
		if not TN+FP == 0:
			Spe = float(TN)/float(TN+FP)
		Acc = float(TP+TN)/float(TP+TN+FP+FN)

		OutFile = JobId+'_p1_vs_p5_PopSize'+str(PopSize)+'_FastMl-MkYlk-f.txt'
		StdMktDf.to_csv(OutFile,sep='\t')

		Comparison = JobId+'_p1_vs_p5_SelThr'+str(SelThr)+'_PopSize'+str(PopSize)
		PredictionDf[Header[0]].append(Comparison)
		PredictionDf[Header[1]].append(TP)
		PredictionDf[Header[2]].append(TN)
		PredictionDf[Header[3]].append(FP)
		PredictionDf[Header[4]].append(FN)
		PredictionDf[Header[5]].append(Sen)
		PredictionDf[Header[6]].append(Spe)
		PredictionDf[Header[7]].append(Acc)

	PredictionDf = pd.DataFrame(PredictionDf,columns=Header)

	return PredictionDf
	

def CreateRocInput():

	## Create ROC input data
	## Input: Mk/Ylk analysis results
	## Output: Input data to plot ROC

	SelThrList = [0]
	MktDfList = sorted(glob.glob('*_FastMl-MkYlk-f.txt'))
	for SelThr in SelThrList:
		for FileName in MktDfList:

			print(FileName)
			RocDf = {};Header = ['ThScore','TP','TN','FP','FN','TPR','FPR']
			for ColName in Header:RocDf.setdefault(ColName,[])
			MktDf = pd.read_csv(FileName,sep='\t',index_col=0)
			MktDf = MktDf.reset_index(drop=True)

			MkScoreList = sorted(MktDf['MkScore'].unique())
			for ThScore in MkScoreList:
				TPR = 'NA';FPR = 'NA'
				TP = len(MktDf[(MktDf['SelCoef'] > SelThr) & (MktDf['MkScore'] >= ThScore)])
				FN = len(MktDf[(MktDf['SelCoef'] > SelThr) & (MktDf['MkScore'] < ThScore)])
				FP = len(MktDf[(MktDf['SelCoef'] <= SelThr) & (MktDf['MkScore'] >= ThScore)])
				TN = len(MktDf[(MktDf['SelCoef'] <= SelThr) & (MktDf['MkScore'] < ThScore)])

				if not TP+FN == 0:TPR = float(TP)/float(TP+FN)
				if not TN+FP == 0:FPR = float(FP)/float(TN+FP)

				RocDf[Header[0]].append(ThScore)
				RocDf[Header[1]].append(TP)
				RocDf[Header[2]].append(TN)
				RocDf[Header[3]].append(FP)
				RocDf[Header[4]].append(FN)
				RocDf[Header[5]].append(TPR)
				RocDf[Header[6]].append(FPR)

			RocDf = pd.DataFrame(RocDf,columns=Header)
			OutFile = FileName.replace('_FastMl-MkYlk-f.txt','_SelThr'+str(SelThr)+'_RocInputDf.txt')
			RocDf.to_csv(OutFile,sep='\t')


if __name__ == '__main__':

	DirectorySetList = [['Scenario-1_None','Scenario-2_P1','Scenario-3_P5','Scenario-4_P1P5']]
	
	BigPredictionDf = pd.DataFrame()
	for DirectorySet in DirectorySetList:

		PopSizeList = ['5','20','100','500']
		for PopSize in PopSizeList:

			p1 = 'Generation_678000_ind_'+PopSize+'_p1_polymorphism_ModDataframe.txt'
			p3 = 'Generation_687000_ind_'+PopSize+'_p3_polymorphism_ModDataframe.txt'
			p4 = 'Generation_760000_ind_'+PopSize+'_p4_polymorphism_ModDataframe.txt'
			p5 = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism_ModDataframe.txt'

			n1 = 'Generation_0_ind_'+PopSize+'_p1_polymorphism_ModDataframe.txt'
			n2 = 'Generation_148000_ind_'+PopSize+'_p1_polymorphism_ModDataframe.txt'
			n3 = 'Generation_491000_ind_'+PopSize+'_p4_polymorphism_ModDataframe.txt'

			PopDic = {'p1':p1, 'p3':p3, 'p4':p4, 'p5':p5,\
					  '1':n1, '2':n2, '3':n3}

			PredictionDf = SummarizeAllSimulation(PopSize,DirectorySet)
			BigPredictionDf = BigPredictionDf.append(PredictionDf)
	
	BigPredictionDf.to_csv('PredictionSummary_MergedScenarios.txt',sep='\t')
	CreateRocInput()
