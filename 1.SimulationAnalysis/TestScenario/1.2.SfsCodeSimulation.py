import glob, os, multiprocessing, re
import argparse
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from scipy import stats

def ExecuteFastMlAsr(IterId,PopSize):

    ## Reconstruct ancestral sequence using FastML
    ## Input: Forward simulation fasta files
    ## Output: '2.2.FastML_PopSize_*/' directories with ancestral sequences

    os.system('mkdir 2.2.FastML_PopSize_'+PopSize)

    TreeFile = 'Simulation_tree.tre'

    RepSeqDic = {}

    P3dir = 'Generation_678000_ind_'+PopSize+'_p3_polymorphism/'
    P2dir = 'Generation_687000_ind_'+PopSize+'_p2_polymorphism/'
    P6dir = 'Generation_760000_ind_'+PopSize+'_p6_polymorphism/'
    P5dir = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism/'

    P3Dic = ParseFastaFile(P3dir+'G0001.fas')
    RepSeqDic['p3'] = P3Dic[list(P3Dic.keys())[0]]
    P2Dic = ParseFastaFile(P2dir+'G0001.fas')
    RepSeqDic['p2'] = P2Dic[list(P2Dic.keys())[0]]
    P6Dic = ParseFastaFile(P6dir+'G0001.fas')
    RepSeqDic['p6'] = P6Dic[list(P6Dic.keys())[0]]
    P5Dic = ParseFastaFile(P5dir+'G0001.fas')
    RepSeqDic['p5'] = P5Dic[list(P5Dic.keys())[0]]

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
                  ' --indelReconstruction BOTH'+\
                  ' --outDir ./ > FastMl.log 2> FastMl.log2'

    #print(Commandline);
    os.system(Commandline)

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

    P3dir = 'Generation_678000_ind_'+PopSize+'_p3_polymorphism/'
    P2dir = 'Generation_687000_ind_'+PopSize+'_p2_polymorphism/'
    P6dir = 'Generation_760000_ind_'+PopSize+'_p6_polymorphism/'
    P5dir = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism/'

    P3Dic = ParseFastaFile(P3dir+'G0001.fas')
    RepSeqDic['p3'] = P3Dic[list(P3Dic.keys())[0]]
    P2Dic = ParseFastaFile(P2dir+'G0001.fas')
    RepSeqDic['p2'] = P2Dic[list(P2Dic.keys())[0]]
    P6Dic = ParseFastaFile(P6dir+'G0001.fas')
    RepSeqDic['p6'] = P6Dic[list(P6Dic.keys())[0]]
    P5Dic = ParseFastaFile(P5dir+'G0001.fas')
    RepSeqDic['p5'] = P5Dic[list(P5Dic.keys())[0]]

    PrankInput = IterId+'.fas'
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
                  ' -F'+\
                  ' -keep '+ \
                  '> 2.1.Prank_PopSize_'+PopSize+'/Prank.log '+ \
                  '2> 2.1.Prank_PopSize_'+PopSize+'/Prank.log2'

    #print(Commandline);
    os.system(Commandline)

    os.system('rm '+PrankInput)


def CreateUnrootedTree(TreeFile):

    UnrootedTree = ''
    fin = open(TreeFile,'r')
    line = fin.readline()
    for i in range(len(line.split('('))):
        if i in [0,2]:
            UnrootedTree += line.split('(')[i]
        else:
            UnrootedTree += '('+line.split('(')[i]
    line = UnrootedTree;UnrootedTree = ''

    for i in range(len(line.split(')'))):

        if i in [range(len(line.split(')')))[-1],range(len(line.split(')')))[-3]]:
            UnrootedTree += line.split(')')[i]
        else:
            UnrootedTree += line.split(')')[i]+')'

    fout = open('UnrootedTree.tre','w')
    fout.write(UnrootedTree);fout.close()


def CreateCtlFile(FileName,PopSize):

    IterId = FileName.split('.')[0]

    fout = open('2.3.Paml_PopSize_'+PopSize+'/'+IterId+'.ctl','w')
    fout.write('	  seqfile = ../'+FileName+'\n')
    fout.write('	 treefile = ../UnrootedTree.tre\n')
    fout.write('	  outfile = '+IterId+'_PAML.txt\n')

    fout.write('		noisy = 0\n')
    fout.write('	  verbose = 0\n')
    fout.write('	  runmode = 0\n')
    fout.write('	  seqtype = 1\n')
    fout.write('	CodonFreq = 2\n')
    fout.write('		clock = 0\n')
    fout.write('		model = 1\n')
    fout.write('	  NSsites = 0\n')
    fout.write('		icode = 0\n')
    fout.write('		Mgene = 0\n')
    fout.write('	fix_kappa = 0\n')
    fout.write('		kappa = 2\n')
    fout.write('	fix_omega = 0\n')
    fout.write('		omega = .4\n')
    fout.write('	fix_alpha = 1\n')
    fout.write('		alpha = 0.\n')
    fout.write('	   Malpha = 0\n')
    fout.write('		ncatG = 8\n')
    fout.write('		getSE = 0\n')
    fout.write(' RateAncestor = 1\n')
    fout.write('	cleandata = 0\n')
    fout.write('	   method = 0\n')
    fout.close()


def ExecutePamlAsr(IterId,PopSize):

    ## Reconstruct ancestral sequence using Paml
    ## Input: Forward simulation fasta files
    ## Output: '2.3.Paml_PopSize_*/' directories with ancestral sequences

    os.system('mkdir 2.3.Paml_PopSize_'+PopSize)

    TreeFile = 'Simulation_tree.tre'
    CreateUnrootedTree(TreeFile)

    RepSeqDic = {}

    P3dir = 'Generation_678000_ind_'+PopSize+'_p3_polymorphism/'
    P2dir = 'Generation_687000_ind_'+PopSize+'_p2_polymorphism/'
    P6dir = 'Generation_760000_ind_'+PopSize+'_p6_polymorphism/'
    P5dir = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism/'

    P3Dic = ParseFastaFile(P3dir+'G0001.fas')
    RepSeqDic['p3'] = P3Dic[list(P3Dic.keys())[0]]
    P2Dic = ParseFastaFile(P2dir+'G0001.fas')
    RepSeqDic['p2'] = P2Dic[list(P2Dic.keys())[0]]
    P6Dic = ParseFastaFile(P6dir+'G0001.fas')
    RepSeqDic['p6'] = P6Dic[list(P6Dic.keys())[0]]
    P5Dic = ParseFastaFile(P5dir+'G0001.fas')
    RepSeqDic['p5'] = P5Dic[list(P5Dic.keys())[0]]

    PamlInput = IterId+'.fas'
    fout = open(PamlInput,'w')
    for s_id in sorted(RepSeqDic.keys()):fout.write('>'+s_id+'\n'+RepSeqDic[s_id]+'\n')
    fout.close()

    CreateCtlFile(PamlInput,PopSize)

    os.chdir('2.3.Paml_PopSize_'+PopSize)

    Commandline = '../../../'+ProgDir+'paml4.9j//bin/codeml '+ \
				  IterId+'.ctl > Paml.log 2> Paml.log2'
    #print(Commandline)
    os.system(Commandline)

    os.chdir('../')


def ReportVariantLoci(GeneId,Sseq,Tseq,SmutDf,TmutDf,PopSize):

    ## Report divergence of sequence between species and polymorphism within population
    ## Input: Subject and target CDS sequences and polymorphisms in dataframe
    ## Output: Dataframe of variants

    MafFilter = float(2)/float(PopSize) # Minimum 2 samples need to have variant
    b_dN = 0;b_pN = 0;dN = 0;dS = 0;pN = 0;pS = 0
    VarLociDf = {};Header = ['IterId','Pos','Sstate','Ref','Alt','SelCoef','GO','VarType']
    for ColName in Header:VarLociDf.setdefault(ColName,[])

    ## Filter for low minor allele frequency
    SpolDf = SmutDf[(SmutDf['Allele_freq'] < 1.0 - MafFilter) & \
                    (SmutDf['Allele_freq'] > 0.0 + MafFilter)].reset_index(drop=True)

    ## Count divergence and polymorphism
    for i in range(len(Sseq)):

        cst = int(i/3)*3;cend = cst+3
        OriCod = Sseq[cst:cend]

        base = i+1
        ## No Divergence ##
        if Sseq[i] == Tseq[i]:

            ## Polymorphism ##
            if len(SpolDf.index) != 0 and base == SpolDf.iloc[0]['MsaPos']:

                OvrlpSpol = SpolDf[SpolDf['MsaPos'] == base]
                for j in OvrlpSpol.index:

                    SelCoef = float(SpolDf.iloc[0]['Sel_coef'])
                    Ref = SpolDf.iloc[0]['Ref'];Alt = SpolDf.iloc[0]['Alt']
                    if Ref == Sseq[i]:
                        VarCod = Sseq[cst:i]+Alt+Sseq[i+1:cend]
                    else:
                        VarCod = Sseq[cst:i]+Ref+Sseq[i+1:cend]

                    VarLociDf[Header[0]].append(GeneId)
                    VarLociDf[Header[1]].append(base)
                    VarLociDf[Header[2]].append(Sseq[i])
                    VarLociDf[Header[3]].append(Ref)
                    VarLociDf[Header[4]].append(Alt)
                    VarLociDf[Header[5]].append(SelCoef)
                    VarLociDf[Header[6]].append(SpolDf.iloc[0]['Emerge_time'])

                    if Seq(OriCod).translate() == Seq(VarCod).translate():
                        pS += 1
                        VarLociDf[Header[7]].append('pS')
                    else:
                        pN += 1
                        VarLociDf[Header[7]].append('pN')
                        if SelCoef > 0:b_pN += 1

                    SpolDf = SpolDf.drop([0]).reset_index(drop=True)

        ## Divergence ##
        elif Sseq[i] != Tseq[i]:
            Smut = SmutDf[(SmutDf['MsaPos'] == base) & (SmutDf['Alt'] == Sseq[i])]

            ## Polymorphism ##
            if len(SpolDf.index) != 0 and base == SpolDf.iloc[0]['MsaPos']:

                OvrlpSpol = SpolDf[SpolDf['MsaPos'] == base]
                States = list(OvrlpSpol['Alt'])+list(OvrlpSpol.iloc[0]['Ref'])
                for j in OvrlpSpol.index:

                    SelCoef = float(SpolDf.iloc[0]['Sel_coef'])
                    Ref = SpolDf.iloc[0]['Ref'];Alt = SpolDf.iloc[0]['Alt']
                    if Ref == Sseq[i]:
                        VarCod = Sseq[cst:i]+Alt+Sseq[i+1:cend]
                    else:
                        VarCod = Sseq[cst:i]+Ref+Sseq[i+1:cend]

                    VarLociDf[Header[0]].append(GeneId)
                    VarLociDf[Header[1]].append(base)
                    VarLociDf[Header[2]].append(Sseq[i])
                    VarLociDf[Header[3]].append(Ref)
                    VarLociDf[Header[4]].append(Alt)
                    VarLociDf[Header[5]].append(SelCoef)
                    VarLociDf[Header[6]].append(SpolDf.iloc[0]['Emerge_time'])

                    if Seq(OriCod).translate() == Seq(VarCod).translate():
                        pS += 1
                        VarLociDf[Header[7]].append('pS')
                    else:
                        pN += 1
                        VarLociDf[Header[7]].append('pN')
                        if SelCoef > 0:b_pN += 1

                    SpolDf = SpolDf.drop([0]).reset_index(drop=True)

                if not Tseq[i] in States:

                    DivCod = Sseq[cst:i]+Tseq[i]+Sseq[i+1:cend]

                    VarLociDf[Header[0]].append(GeneId)
                    VarLociDf[Header[1]].append(base)
                    VarLociDf[Header[2]].append(Sseq[i])
                    VarLociDf[Header[3]].append(Sseq[i])
                    VarLociDf[Header[4]].append(Tseq[i])

                    SelCoef = 0;GO = 0
                    if len(Smut.index) != 0:
                        SelCoef = float(Smut['Sel_coef'].values[0])
                        GO = int(Smut['Emerge_time'].values[0])

                    VarLociDf[Header[5]].append(SelCoef)
                    VarLociDf[Header[6]].append(GO)

                    if Seq(OriCod).translate() == Seq(DivCod).translate():
                        dS += 1
                        VarLociDf[Header[7]].append('dS')
                    else:
                        dN += 1
                        VarLociDf[Header[7]].append('dN')
                        if SelCoef > 0:b_dN += 1

            ## No polymorphism ##
            else:

                DivCod = Sseq[cst:i]+Tseq[i]+Sseq[i+1:cend]

                VarLociDf[Header[0]].append(GeneId)
                VarLociDf[Header[1]].append(base)
                VarLociDf[Header[2]].append(Sseq[i])
                VarLociDf[Header[3]].append(Sseq[i])
                VarLociDf[Header[4]].append(Tseq[i])

                SelCoef = 0;GO = 0;
                if len(Smut.index) != 0:
                    SelCoef = float(Smut['Sel_coef'].values[0])
                    GO = int(Smut['Emerge_time'].values[0])

                VarLociDf[Header[5]].append(SelCoef)
                VarLociDf[Header[6]].append(GO)
                if Seq(OriCod).translate() == Seq(DivCod).translate():
                    dS += 1
                    VarLociDf[Header[7]].append('dS')
                else:
                    dN += 1
                    VarLociDf[Header[7]].append('dN')
                    if SelCoef > 0:b_dN += 1

    VarLociDf = pd.DataFrame(VarLociDf,columns=Header)

    return VarLociDf, dN,dS,pN,pS,b_dN,b_pN


def ParseFastaFile(fname):

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

        MutDf = pd.read_csv(f,sep='\t',index_col=0)
        MutDf = MutDf.sort_values(by=['Position']).reset_index(drop=True)
        for i in range(len(MutDf.index)):

            Pos = MutDf.iloc[i]['Position'] - 1

            while Pos >= CurEnd and len(GeneCoDf.index) != 0:

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


def PrankMknMK(IterId,Subject,Target, PopSize):

    ## Count divergences (dN, dS) and polymorphisms (pN, pS) and calculate MK score
    ## Input: Variant calling files and Fasta file containing aligned sequences
    ## Output: Variant report dataframe and Mk/nMK dataframe

    IdDic = {'p3':'p3', 'p2':'p2', 'p6':'p6', 'p5':'p5', '#3#':'p0', '#2#':'p1', '#1#':'p4'}

    p3 = 'Generation_678000_ind_'+PopSize+'_p3_polymorphism.vcf'
    p2 = 'Generation_687000_ind_'+PopSize+'_p2_polymorphism.vcf'
    p6 = 'Generation_760000_ind_'+PopSize+'_p6_polymorphism.vcf'
    p5 = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism.vcf'

    p0 = 'Generation_0_ind_'+PopSize+'_p0_polymorphism.vcf'
    p1 = 'Generation_148000_ind_'+PopSize+'_p1_polymorphism.vcf'
    p4 = 'Generation_491000_ind_'+PopSize+'_p4_polymorphism.vcf'

    PopDic = {'p3':p3, 'p2':p2, 'p6':p6, 'p5':p5, '#3#':p0, '#2#':p1, '#1#':p4}

    SubjectVcf = PopDic[Subject];TargetVcf = PopDic[Target]
    SmutDf = pd.read_csv(SubjectVcf.replace('.vcf','_ModDataframe.txt'),sep='\t',index_col=0)
    TmutDf = pd.read_csv(TargetVcf.replace('.vcf','_ModDataframe.txt'),sep='\t',index_col=0)

    FileName = '2.1.Prank_PopSize_'+PopSize+'/'+IterId+'.anc.fas'
    SeqDic = ParseFastaFile(FileName)
    Sseq = SeqDic[Subject]
    Tseq = SeqDic[Target]

    MktDf = {};MktHeader = ['IterId','dN','dS','pN','pS']
    for ColName in MktHeader:MktDf.setdefault(ColName,[])

    GeneVarDf,dN,dS,pN,pS,bdN,bpN = ReportVariantLoci(IterId,Sseq,Tseq,SmutDf,TmutDf,PopSize)

    MktDf[MktHeader[0]].append(IterId)
    MktDf[MktHeader[1]].append(dN)
    MktDf[MktHeader[2]].append(dS)
    MktDf[MktHeader[3]].append(pN)
    MktDf[MktHeader[4]].append(pS)

    Subject = IdDic[Subject];Target = IdDic[Target]
    MktDf = pd.DataFrame(MktDf,columns=MktHeader)
    MktDf.to_csv(Subject+'_vs_'+Target+'_Popsize'+PopSize+'_Prank-MknMK.txt',sep='\t')
    GeneVarDf.to_csv(Subject+'_vs_'+Target+'_Popsize'+PopSize+'_Prank-VarRept.txt',sep='\t')


def FastMlMknMK(IterId,Subject,Target, PopSize):

    ## Count divergences (dN, dS) and polymorphisms (pN, pS) and calculate MK score
    ## Input: Variant calling files and Fasta file containing aligned sequences
    ## Output: Variant report dataframe and Mk/nMK dataframe

    IdDic = {'p3':'p3', 'p2':'p2', 'p6':'p6', 'p5':'p5', '1':'p0', '2':'p1', '3':'p4'}

    p3 = 'Generation_678000_ind_'+PopSize+'_p3_polymorphism.vcf'
    p2 = 'Generation_687000_ind_'+PopSize+'_p2_polymorphism.vcf'
    p6 = 'Generation_760000_ind_'+PopSize+'_p6_polymorphism.vcf'
    p5 = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism.vcf'

    p0 = 'Generation_0_ind_'+PopSize+'_p0_polymorphism.vcf'
    p1 = 'Generation_148000_ind_'+PopSize+'_p1_polymorphism.vcf'
    p4 = 'Generation_491000_ind_'+PopSize+'_p4_polymorphism.vcf'

    PopDic = {'p3':p3, 'p2':p2, 'p6':p6, 'p5':p5, '1':p0, '2':p1, '3':p4}

    SubjectVcf = PopDic[Subject];TargetVcf = PopDic[Target]
    SmutDf = pd.read_csv(SubjectVcf.replace('.vcf','_ModDataframe.txt'),sep='\t',index_col=0)
    TmutDf = pd.read_csv(TargetVcf.replace('.vcf','_ModDataframe.txt'),sep='\t',index_col=0)

    FileName = '2.2.FastML_PopSize_'+PopSize+'/FilesForJalView/seq.joint_JalView.FASTA.aln'
    SeqDic = ParseFastaFile(FileName)
    Sseq = SeqDic[Subject]
    Tseq = SeqDic[Target]

    MktDf = {};MktHeader = ['IterId','dN','dS','pN','pS']
    for ColName in MktHeader:MktDf.setdefault(ColName,[])

    GeneVarDf,dN,dS,pN,pS,bdN,bpN = ReportVariantLoci(IterId,Sseq,Tseq,SmutDf,TmutDf,PopSize)

    MktDf[MktHeader[0]].append(IterId)
    MktDf[MktHeader[1]].append(dN)
    MktDf[MktHeader[2]].append(dS)
    MktDf[MktHeader[3]].append(pN)
    MktDf[MktHeader[4]].append(pS)

    Subject = IdDic[Subject];Target = IdDic[Target]
    MktDf = pd.DataFrame(MktDf,columns=MktHeader)
    MktDf.to_csv(Subject+'_vs_'+Target+'_Popsize'+PopSize+'_FastMl-MknMK.txt',sep='\t')
    GeneVarDf.to_csv(Subject+'_vs_'+Target+'_Popsize'+PopSize+'_Fastml-VarRept.txt',sep='\t')


def ParsePamlAsr(FileName):

    SeqDic = {};flag = 0;cnt = 0
    fin = open(FileName,'r')
    for line in fin:
        if 'List of extant and reconstructed sequences' in line:
            flag = 1
        if flag == 1:

            cnt += 1
            if cnt == 3:Length = int(line.strip().split()[0])
            if cnt >= 5 and cnt < 5+Length:

                SeqId = re.split(r'\s{2,}',line.strip())[0]
                Sequence = re.split(r'\s{2,}',line.strip())[1].replace(' ','')
                SeqDic[SeqId] = Sequence

    fin.close()

    return SeqDic


def PamlMknMK(IterId,Subject,Target, PopSize):
	
    ## Count divergences (dN, dS) and polymorphisms (pN, pS) and calculate MK score
    ## Input: Variant calling files and Fasta file containing aligned sequences
    ## Output: Variant report dataframe and Mk/nMK dataframe

    IdDic = {'p3':'p3', 'p2':'p2', 'p6':'p6', 'p5':'p5', 'node #5':'p1', 'node #6':'p4'}

    p3 = 'Generation_678000_ind_'+PopSize+'_p3_polymorphism.vcf'
    p2 = 'Generation_687000_ind_'+PopSize+'_p2_polymorphism.vcf'
    p6 = 'Generation_760000_ind_'+PopSize+'_p6_polymorphism.vcf'
    p5 = 'Generation_665000_ind_'+PopSize+'_p5_polymorphism.vcf'

    p0 = 'Generation_0_ind_'+PopSize+'_p0_polymorphism.vcf'
    p1 = 'Generation_148000_ind_'+PopSize+'_p1_polymorphism.vcf'
    p4 = 'Generation_491000_ind_'+PopSize+'_p4_polymorphism.vcf'

    PopDic = {'p3':p3, 'p2':p2, 'p6':p6, 'p5':p5, 'node #5':p1, 'node #6':p4}

    SubjectVcf = PopDic[Subject];TargetVcf = PopDic[Target]
    SmutDf = pd.read_csv(SubjectVcf.replace('.vcf','_ModDataframe.txt'),sep='\t',index_col=0)
    TmutDf = pd.read_csv(TargetVcf.replace('.vcf','_ModDataframe.txt'),sep='\t',index_col=0)

    FileName = '2.3.Paml_PopSize_'+PopSize+'/rst'
    SeqDic = ParsePamlAsr(FileName)
    Sseq = SeqDic[Subject]
    Tseq = SeqDic[Target]

    MktDf = {};MktHeader = ['IterId','dN','dS','pN','pS']
    for ColName in MktHeader:MktDf.setdefault(ColName,[])

    GeneVarDf,dN,dS,pN,pS,bdN,bpN = ReportVariantLoci(IterId,Sseq,Tseq,SmutDf,TmutDf,PopSize)

    MktDf[MktHeader[0]].append(IterId)
    MktDf[MktHeader[1]].append(dN)
    MktDf[MktHeader[2]].append(dS)
    MktDf[MktHeader[3]].append(pN)
    MktDf[MktHeader[4]].append(pS)

    Subject = IdDic[Subject];Target = IdDic[Target]
    MktDf = pd.DataFrame(MktDf,columns=MktHeader)
    MktDf.to_csv(Subject+'_vs_'+Target+'_Popsize'+PopSize+'_Paml-MknMK.txt',sep='\t')
    GeneVarDf.to_csv(Subject+'_vs_'+Target+'_Popsize'+PopSize+'_Paml-VarRept.txt',sep='\t')


def ExecuteMultipleSimulation(Directory,SimNum,ThrNum):

    ## Perform Simulation in parallel
    ## Simulation options can be found in 'TemplateData/' Directory
    ## Input: Simulation Directory containing 'TemplateData/' Directory
    ## Output: Simulation data representing population of individuals their genes

    os.chdir(Directory);IterDirList = []

    for i in range(SimNum):IterDirList.append(Directory+'/Iter'+str(i+1).zfill(4))

    p = multiprocessing.Pool(ThrNum)
    p.map(ExecuteSfsCode,IterDirList)

    os.chdir('../')


def ExecuteSfsCode(Directory):

    ## Perform Simulation using SFS_CODE
    ## Input: Copies of 'TemplateData/' directory containing simulation eidos files
    ## Output: Simulation data and Mk/nMK analysis results

    IterId = Directory.split('/')[-1]
    if not IterId in os.listdir('./'):
    	os.system('cp -R TemplateData '+IterId)
    os.chdir(IterId)
	
    RunId = Directory.split('/')[-2]+' '+IterId
    systemstr = 'python ExecuteSfsCode.py'
    print(RunId+' SFS_Code Simulation');os.system(systemstr)

    for PopSize in PopSizeList:

        systemstr = 'python 001.Output_fasta_file_ver4.py '+PopSize
        print(RunId+' Create Fasta File, PopSize = '+PopSize);os.system(systemstr)

        ModifyMutDf(PopSize)

        if AncMet == 'prank':

            ExecutePrankAsr(Directory.split('/')[-1],PopSize)
            PrankMknMK(Directory.split('/')[-1],'p3','#2#',PopSize)

        elif AncMet == 'fastml':

            ExecuteFastMlAsr(Directory.split('/')[-1],PopSize)
            FastMlMknMK(Directory.split('/')[-1],'p3','2',PopSize)

        elif AncMed == 'paml':

            ExecutePamlAsr(Directory.split('/')[-1],PopSize)
            PamlMknMK(Directory.split('/')[-1],'p3','node #5',PopSize)


        ## Execute standard MK
        if (args.std):

            FastMlMknMK(Directory.split('/')[-1],'p3','p5',PopSize) #

    os.chdir('../')


def MakePrediction(MknMKDf,MkThresholds):

    PosScens = ['Scenario-1_+_+','Scenario-2_+_0','Scenario-3_+_-']
    NeuScens = ['Scenario-4_0_+','Scenario-5_0_0','Scenario-6_0_-']
    NegScens = ['Scenario-7_-_+','Scenario-8_-_0','Scenario-9_-_-']

    Divergence = 29600
    SelCoefList = [];PosSelList = []
    PosPos = [0,0,0,0,0];PosNeu = [0,0,0,0,0];PosNeg = [0,0,0,0,0]
    NeuPos = [0,0,0,0,0];NeuNeu = [0,0,0,0,0];NeuNeg = [0,0,0,0,0]
    NegPos = [0,0,0,0,0];NegNeu = [0,0,0,0,0];NegNeg = [0,0,0,0,0]

    p3 = 'Generation_678000_ind_100_p3_polymorphism_ModDataframe.txt'
    p2 = 'Generation_687000_ind_100_p2_polymorphism_ModDataframe.txt'
    p6 = 'Generation_760000_ind_100_p6_polymorphism_ModDataframe.txt'
    p5 = 'Generation_665000_ind_100_p5_polymorphism_ModDataframe.txt'

    p0 = 'Generation_0_ind_100_p0_polymorphism_ModDataframe.txt'
    p1 = 'Generation_148000_ind_100_p1_polymorphism_ModDataframe.txt'
    p4 = 'Generation_491000_ind_100_p4_polymorphism_ModDataframe.txt'

    PopDic = {'p3':p3, 'p2':p2, 'p6':p6, 'p5':p5, 'p0':p0, 'p1':p1, 'p4':p4}

    for i in MknMKDf.index:

        IterId = MknMKDf.loc[i,'IterId']
        MktScore = [MknMKDf.loc[i,'MkScore_1'],MknMKDf.loc[i,'MkScore_2'],MknMKDf.loc[i,'MkScore_3'],MknMKDf.loc[i,'MkScore_4'],MknMKDf.loc[i,'MkScore_5']]
        Scenario = MknMKDf.loc[i,'Scenario']
        VarDfFile = PopDic['p3']
        VarDf = pd.read_csv(IterId+'/'+VarDfFile,sep='\t',index_col=0)
        RedVarDf = VarDf[(VarDf['Emerge_time'] >= Divergence) & (VarDf['Allele_freq'] ==1.0)]

        for i in range(5):
            if Scenario in PosScens:
                PosSel = '+'
                if MktScore[i] > MkThresholds[0][i]:
                    PosPos[i] += 1
                elif MktScore[i] < MkThresholds[1][i]:
                    PosNeg[i] += 1
                else:
                    PosNeu[i] += 1
            elif Scenario in NegScens:
                PosSel = '-'
                if MktScore[i] > MkThresholds[0][i]:
                    NegPos[i] += 1
                elif MktScore[i] < MkThresholds[1][i]:
                    NegNeg[i] += 1
                else:
                    NegNeu[i] += 1
            else:
                PosSel = '0'
                if MktScore[i] > MkThresholds[0][i]:
                    NeuPos[i] += 1
                elif MktScore[i] < MkThresholds[1][i]:
                    NeuNeg[i] += 1
                else:
                    NeuNeu[i] += 1

        SelCoefList.append(sum(RedVarDf['Sel_coef']))
        PosSelList.append(PosSel)

    Results = [PosPos,PosNeu,PosNeg,NeuPos,NeuNeu,NeuNeg,NegPos,NegNeu,NegNeg]

    return Results,SelCoefList,PosSelList


def SummarizeAllSimulation(PopSize,Directory):

    ## Summarize results
    ## Input: Simulation data and Mk/nMK analysis results
    ## Output: Performance of Mk/nMK analyses

    # Calculate expected [dN dS pN pS] & Mk score threshold
    IterList = sorted(glob.glob('Scenario-5_0_0/Iter*'))
    NeutralDf = pd.DataFrame()
    for IterId in IterList:

        if AncMet == 'prank':
            MktFile = IterId+'/p3_vs_p1_Popsize'+PopSize+'_Prank-MknMK.txt'
        elif AncMet == 'fastml':
            MktFile = IterId+'/p3_vs_p1_Popsize'+PopSize+'_FastMl-MknMK.txt'
        elif AncMet == 'paml':
            MktFile = IterId+'/p3_vs_p1_Popsize'+PopSize+'_Paml-MknMK.txt'

        SubDf = pd.read_csv(MktFile,sep='\t',index_col=0)
        NeutralDf = NeutralDf.append(SubDf)

    AvgdN = NeutralDf['dN'].mean();AvgdS = NeutralDf['dS'].mean()
    AvgpN = NeutralDf['pN'].mean();AvgpS = NeutralDf['pS'].mean()

    QuantileThr = 0.1

    # Ver 1
    NeutralMkScores_1 = ((NeutralDf['dN']+1)/(NeutralDf['dS']+1))/((NeutralDf['pN']+1)/(NeutralDf['pS']+1))

    # Ver 2
    NeutralMkScores_2 = ((NeutralDf['dN']+AvgdN)/(NeutralDf['dS']+AvgdS))/((NeutralDf['pN']+AvgpN)/(NeutralDf['pS']+AvgpS))

    # Ver 3
    NeutralMkScores_3 = np.log(np.exp(((NeutralDf['dN']+AvgdN)/AvgdN)-((NeutralDf['dS']+AvgdS)/AvgdS))/np.exp(((NeutralDf['pN']+AvgpN)/AvgpN)-((NeutralDf['pS']+AvgpS)/AvgpS)))

    # Ver 4
    NeutralMkScores_4 = np.log(np.exp(((NeutralDf['dN']+AvgdN)/AvgdN)**2-((NeutralDf['dS']+AvgdS)/AvgdS)**2)/np.exp(((NeutralDf['pN']+AvgpN)/AvgpN)**2-((NeutralDf['pS']+AvgpS)/AvgpS)**2))

    # Ver 5
    AvgDiv = sum([AvgdN,AvgdS]);AvgPol = sum([AvgpN,AvgpS])
    NeutralMkScores_5 = np.log(np.exp(((NeutralDf['dN']+AvgdN)/AvgDiv)**2-((NeutralDf['dS']+AvgdS)/AvgDiv)**2)/np.exp(((NeutralDf['pN']+AvgpN)/AvgPol)**2-((NeutralDf['pS']+AvgpS)/AvgPol)**2))
    
    MkThreshold_Pos_1 = NeutralMkScores_1.quantile(1-QuantileThr)
    MkThreshold_Neg_1 = NeutralMkScores_1.quantile(QuantileThr)
    MkThreshold_Pos_2 = NeutralMkScores_2.quantile(1-QuantileThr)
    MkThreshold_Neg_2 = NeutralMkScores_2.quantile(QuantileThr)
    MkThreshold_Pos_3 = NeutralMkScores_3.quantile(1-QuantileThr)
    MkThreshold_Neg_3 = NeutralMkScores_3.quantile(QuantileThr)
    MkThreshold_Pos_4 = NeutralMkScores_4.quantile(1-QuantileThr)
    MkThreshold_Neg_4 = NeutralMkScores_4.quantile(QuantileThr)
    MkThreshold_Pos_5 = NeutralMkScores_5.quantile(1-QuantileThr)
    MkThreshold_Neg_5 = NeutralMkScores_5.quantile(QuantileThr)

    MkThresholds = [[MkThreshold_Pos_1,MkThreshold_Pos_2,MkThreshold_Pos_3,MkThreshold_Pos_4,MkThreshold_Pos_5], \
                    [MkThreshold_Neg_1,MkThreshold_Neg_2,MkThreshold_Neg_3,MkThreshold_Neg_4,MkThreshold_Neg_5]]

    os.chdir(Directory)

    PredictionDf = {}
    Header = ['Method','Selection','PopSize', \
            'PosPos','PosNeu','PosNeg','NeuPos','NeuNeu','NeuNeg','NegPos','NegNeu','NegNeg', \
            'PosSen','PosSpe','PosAcc','NegSen','NegSpe','NegAcc']
    for ColName in Header:PredictionDf.setdefault(ColName,[])

    IterList = sorted(glob.glob('Iter*'))

    ## Novel MK ##
    MKnMKDf = pd.DataFrame()
    for Element in IterList:

        IterId = Element.split('/')[-1]

        if AncMet == 'prank':
            MktFile = IterId+'/p3_vs_p1_Popsize'+PopSize+'_Prank-MknMK.txt'
        elif AncMet == 'fastml':
            MktFile = IterId+'/p3_vs_p1_Popsize'+PopSize+'_FastMl-MknMK.txt'
        elif AncMet == 'paml':
            MktFile = IterId+'/p3_vs_p1_Popsize'+PopSize+'_Paml-MknMK.txt'

        SubDf = pd.read_csv(MktFile,sep='\t',index_col=0)

        MkScore_1 = ((SubDf.iloc[0]['dN']+1)/(SubDf.iloc[0]['dS']+1))/((SubDf.iloc[0]['pN']+1)/(SubDf.iloc[0]['pS']+1))
        MkScore_2 = ((SubDf.iloc[0]['dN']+AvgdN)/(SubDf.iloc[0]['dS']+AvgdS))/((SubDf.iloc[0]['pN']+AvgpN)/(SubDf.iloc[0]['pS']+AvgpS))
        MkScore_3 = np.log(np.exp(((SubDf.iloc[0]['dN']+AvgdN)/AvgdN)-((SubDf.iloc[0]['dS']+AvgdS)/AvgdS))/np.exp(((SubDf.iloc[0]['pN']+AvgpN)/AvgpN)-((SubDf.iloc[0]['pS']+AvgpS)/AvgpS)))
        MkScore_4 = np.log(np.exp(((SubDf.iloc[0]['dN']+AvgdN)/AvgdN)**2-((SubDf.iloc[0]['dS']+AvgdS)/AvgdS)**2)/np.exp(((SubDf.iloc[0]['pN']+AvgpN)/AvgpN)**2-((SubDf.iloc[0]['pS']+AvgpS)/AvgpS)**2))
        AvgDiv = sum([AvgdN,AvgdS]);AvgPol = sum([AvgpN,AvgpS])
        MkScore_5 = np.log(np.exp(((SubDf.iloc[0]['dN']+AvgdN)/AvgDiv)**2-((SubDf.iloc[0]['dS']+AvgdS)/AvgDiv)**2)/np.exp(((SubDf.iloc[0]['pN']+AvgpN)/AvgPol)**2-((SubDf.iloc[0]['pS']+AvgpS)/AvgPol)**2))

        SubDf['Scenario'] = Directory.replace('/','')
        SubDf['MkScore_1'] = MkScore_1
        SubDf['MkScore_2'] = MkScore_2
        SubDf['MkScore_3'] = MkScore_3
        SubDf['MkScore_4'] = MkScore_4
        SubDf['MkScore_5'] = MkScore_5

        MKnMKDf = MKnMKDf.append(SubDf)

    MKnMKDf = MKnMKDf.reset_index(drop=True)

    Results,SelCoefList,PosSelList = MakePrediction(MKnMKDf,MkThresholds)
    PosPos,PosNeu,PosNeg,NeuPos,NeuNeu,NeuNeg,NegPos,NegNeu,NegNeg = Results # Prediction results

    # Sensitivity, Specificity and Accuracy for positively selected regions
    PosSens = '';PosSpes = '';PosAccs = ''
    for i in range(len(PosPos)):
        TP = PosPos[i];FN = PosNeu[i]+PosNeg[i]
        FP = NeuPos[i]+NegPos[i];TN = NeuNeu[i]+NeuNeg[i]+NegNeu[i]+NegNeg[i]
        PosSen = 'NA';PosSpe = 'NA';PosAcc = 'NA'
        if not TP+FN == 0:
            PosSen = float(TP)/float(TP+FN)
            PosSens += str(PosSen)+';'
        if not TN+FP == 0:
            PosSpe = float(TN)/float(TN+FP)
            PosSpes += str(PosSpe)+';'
        PosAcc = float(TP+TN)/float(TP+TN+FP+FN)
        PosAccs += str(PosAcc)+';'

    # Sensitivity, Specificity and Accuracy for negatively selected regions
    NegSens = '';NegSpes = '';NegAccs = ''

    for i in range(len(NegNeg)):
        TP = NegNeg[i];FN = NegPos[i]+NegNeu[i]
        FP = PosNeg[i]+NeuNeg[i];TN = PosPos[i]+PosNeu[i]+NeuPos[i]+NeuNeu[i]
        NegSen = 'NA';NegSpe = 'NA';NegAcc = 'NA'
        if not TP+FN == 0:
            NegSen = float(TP)/float(TP+FN)
            NegSens += str(NegSen)+';'
        if not TN+FP == 0:
            NegSpe = float(TN)/float(TN+FP)
            NegSpes += str(NegSpe)+';'
        NegAcc = float(TP+TN)/float(TP+TN+FP+FN)
        NegAccs += str(NegAcc)+';'

    if AncMet == 'prank':
        OutFile = 'p3_vs_p1_PopSize'+str(PopSize)+'_Prank-MknMK-f.txt'
    elif AncMet == 'fastml':
        OutFile = 'p3_vs_p1_PopSize'+str(PopSize)+'_FastMl-MknMK-f.txt'
    elif AncMet == 'paml':
        OutFile = 'p3_vs_p1_PopSize'+str(PopSize)+'_Paml-MknMK-f.txt'

    MKnMKDf.to_csv(OutFile,sep='\t')

    Elements = ['nMK-'+AncMet,Directory,PopSize, \
                PosPos,PosNeu,PosNeg,NeuPos,NeuNeu,NeuNeg,NegPos,NegNeu,NegNeg, \
                PosSens,PosSpes,PosAccs,NegSens,NegSpes,NegAccs]

    for i in range(len(Elements)):PredictionDf[Header[i]].append(Elements[i])

    ## Standard MK ##
    if (args.std):

        # Calculate expected [dN dS pN pS] & Mk score threshold
        IterList = sorted(glob.glob('../Scenario-5_0_0/Iter*'))
        NeutralDf = pd.DataFrame()
        for IterId in IterList:

            MktFile = IterId+'/p3_vs_p5_Popsize'+PopSize+'_FastMl-MknMK.txt'
            SubDf = pd.read_csv(MktFile,sep='\t',index_col=0)
            NeutralDf = NeutralDf.append(SubDf)

        NeutralMkScores_1 = ((NeutralDf['dN']+1)/(NeutralDf['dS']+1))/((NeutralDf['pN']+1)/(NeutralDf['pS']+1))
        NeutralMkScores_2 = ((NeutralDf['dN']+AvgdN)/(NeutralDf['dS']+AvgdS))/((NeutralDf['pN']+AvgpN)/(NeutralDf['pS']+AvgpS))
        NeutralMkScores_3 = np.log(np.exp(((NeutralDf['dN']+AvgdN)/AvgdN)-((NeutralDf['dS']+AvgdS)/AvgdS))/np.exp(((NeutralDf['pN']+AvgpN)/AvgpN)-((NeutralDf['pS']+AvgpS)/AvgpS)))
        NeutralMkScores_4 = np.log(np.exp(((NeutralDf['dN']+AvgdN)/AvgdN)**2-((NeutralDf['dS']+AvgdS)/AvgdS)**2)/np.exp(((NeutralDf['pN']+AvgpN)/AvgpN)**2-((NeutralDf['pS']+AvgpS)/AvgpS)**2))
        AvgDiv = sum([AvgdN,AvgdS]);AvgPol = sum([AvgpN,AvgpS])
        NeutralMkScores_5 = np.log(np.exp(((NeutralDf['dN']+AvgdN)/AvgDiv)**2-((NeutralDf['dS']+AvgdS)/AvgDiv)**2)/np.exp(((NeutralDf['pN']+AvgpN)/AvgPol)**2-((NeutralDf['pS']+AvgpS)/AvgPol)**2))

        MkThreshold_Pos_1 = NeutralMkScores_1.quantile(1-QuantileThr)
        MkThreshold_Neg_1 = NeutralMkScores_1.quantile(QuantileThr)
        MkThreshold_Pos_2 = NeutralMkScores_2.quantile(1-QuantileThr)
        MkThreshold_Neg_2 = NeutralMkScores_2.quantile(QuantileThr)
        MkThreshold_Pos_3 = NeutralMkScores_3.quantile(1-QuantileThr)
        MkThreshold_Neg_3 = NeutralMkScores_3.quantile(QuantileThr)
        MkThreshold_Pos_4 = NeutralMkScores_4.quantile(1-QuantileThr)
        MkThreshold_Neg_4 = NeutralMkScores_4.quantile(QuantileThr)
        MkThreshold_Pos_5 = NeutralMkScores_5.quantile(1-QuantileThr)
        MkThreshold_Neg_5 = NeutralMkScores_5.quantile(QuantileThr)

        MkThresholds = [[MkThreshold_Pos_1,MkThreshold_Pos_2,MkThreshold_Pos_3,MkThreshold_Pos_4,MkThreshold_Pos_5], \
                       [MkThreshold_Neg_1,MkThreshold_Neg_2,MkThreshold_Neg_3,MkThreshold_Neg_4,MkThreshold_Neg_5]]

        IterList = sorted(glob.glob('Iter*'))

        StdMktDf = pd.DataFrame()
        for IterId in IterList:
		
            MktFile = IterId+'/p3_vs_p5_Popsize'+PopSize+'_FastMl-MknMK.txt'
            SubDf = pd.read_csv(MktFile,sep='\t',index_col=0)

            MkScore_1 = ((SubDf.iloc[0]['dN']+1)/(SubDf.iloc[0]['dS']+1))/((SubDf.iloc[0]['pN']+1)/(SubDf.iloc[0]['pS']+1))
            MkScore_2 = ((SubDf.iloc[0]['dN']+AvgdN)/(SubDf.iloc[0]['dS']+AvgdS))/((SubDf.iloc[0]['pN']+AvgpN)/(SubDf.iloc[0]['pS']+AvgpS))
            MkScore_3 = np.log(np.exp(((SubDf.iloc[0]['dN']+AvgdN)/AvgdN)-((SubDf.iloc[0]['dS']+AvgdS)/AvgdS))/np.exp(((SubDf.iloc[0]['pN']+AvgpN)/AvgpN)-((SubDf.iloc[0]['pS']+AvgpS)/AvgpS)))
            MkScore_4 = np.log(np.exp(((SubDf.iloc[0]['dN']+AvgdN)/AvgdN)**2-((SubDf.iloc[0]['dS']+AvgdS)/AvgdS)**2)/np.exp(((SubDf.iloc[0]['pN']+AvgpN)/AvgpN)**2-((SubDf.iloc[0]['pS']+AvgpS)/AvgpS)**2))
            AvgDiv = sum([AvgdN,AvgdS]);AvgPol = sum([AvgpN,AvgpS])
            MkScore_5 = np.log(np.exp(((SubDf.iloc[0]['dN']+AvgdN)/AvgDiv)**2-((SubDf.iloc[0]['dS']+AvgdS)/AvgDiv)**2)/np.exp(((SubDf.iloc[0]['pN']+AvgpN)/AvgPol)**2-((SubDf.iloc[0]['pS']+AvgpS)/AvgPol)**2))

            SubDf['Scenario'] = Directory.replace('/','')
            SubDf['MkScore_1'] = MkScore_1
            SubDf['MkScore_2'] = MkScore_2
            SubDf['MkScore_3'] = MkScore_3
            SubDf['MkScore_4'] = MkScore_4
            SubDf['MkScore_5'] = MkScore_5

            StdMktDf = StdMktDf.append(SubDf)
	
        StdMktDf = StdMktDf.reset_index(drop=True)
        Results,SelCoefList,PosSelList = MakePrediction(StdMktDf,MkThresholds)
        PosPos,PosNeu,PosNeg,NeuPos,NeuNeu,NeuNeg,NegPos,NegNeu,NegNeg = Results # Prediction results


        # Sensitivity, Specificity and Accuracy for positively selected regions
        PosSens = '';PosSpes = '';PosAccs = ''
        for i in range(len(PosPos)):
            TP = PosPos[i];FN = PosNeu[i]+PosNeg[i]
            FP = NeuPos[i]+NegPos[i];TN = NeuNeu[i]+NeuNeg[i]+NegNeu[i]+NegNeg[i]
            PosSen = 'NA';PosSpe = 'NA';PosAcc = 'NA'
            if not TP+FN == 0:
                PosSen = float(TP)/float(TP+FN)
                PosSens += str(PosSen)+';'
            if not TN+FP == 0:
                PosSpe = float(TN)/float(TN+FP)
                PosSpes += str(PosSpe)+';'
            PosAcc = float(TP+TN)/float(TP+TN+FP+FN)
            PosAccs += str(PosAcc)+';'

        # Sensitivity, Specificity and Accuracy for negatively selected regions
        NegSens = '';NegSpes = '';NegAccs = ''
        for i in range(len(NegNeg)):
            TP = NegNeg[i];FN = NegPos[i]+NegNeu[i]
            FP = PosNeg[i]+NeuNeg[i];TN = PosPos[i]+PosNeu[i]+NeuPos[i]+NeuNeu[i]
            NegSen = 'NA';NegSpe = 'NA';NegAcc = 'NA'
            if not TP+FN == 0:
                NegSen = float(TP)/float(TP+FN)
                NegSens += str(NegSen)+';'
            if not TN+FP == 0:
                NegSpe = float(TN)/float(TN+FP)
                NegSpes += str(NegSpe)+';'
            NegAcc = float(TP+TN)/float(TP+TN+FP+FN)
            NegAccs += str(NegAcc)+';'

        OutFile = 'p3_vs_p5_PopSize'+str(PopSize)+'_MknMK-f.txt'
        StdMktDf.to_csv(OutFile,sep='\t')

        Elements = ['Std. MK',Directory,PopSize, \
                    PosPos,PosNeu,PosNeg,NeuPos,NeuNeu,NeuNeg,NegPos,NegNeu,NegNeg, \
                    PosSens,PosSpes,PosAccs,NegSens,NegSpes,NegAccs]
        for i in range(len(Elements)):PredictionDf[Header[i]].append(Elements[i])

    PredictionDf = pd.DataFrame(PredictionDf,columns=Header)

    os.chdir('../')

    return PredictionDf


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--ancmet', choices=['prank','fastml','paml'], required=True,
            metavar='[Program name]',
            help='Name of the ancestral sequence reconstruction method used')
    parser.add_argument('--std', action='store_true')
	
    args = parser.parse_args()

    AncMet = args.ancmet

    MainDir = './'
    ProgDir = '../Program/'
    SimNum = 1000
    ThrNum = 50

    BigPredictionDf = pd.DataFrame()
    ScenarioList = ['Scenario-1_+_+','Scenario-2_+_0','Scenario-3_+_-', \
                    'Scenario-4_0_+','Scenario-5_0_0','Scenario-6_0_-', \
                    'Scenario-7_-_+','Scenario-8_-_0','Scenario-9_-_-']
    #ScenarioList = ['Scenario-4_0_+','Scenario-6_0_-', \
    # 				'Scenario-7_-_+','Scenario-8_-_0']
    # ScenarioList = ['Scenario-9_-_-']
    PopSizeList = ['5','20','100']
	
    #for Directory in ScenarioList:

    #	ExecuteMultipleSimulation(Directory,SimNum,ThrNum)

	
    for Directory in ScenarioList:

        for PopSize in PopSizeList:

            PredictionDf = SummarizeAllSimulation(PopSize,Directory)
            BigPredictionDf = BigPredictionDf.append(PredictionDf)

    BigPredictionDf = BigPredictionDf.reset_index(drop=True)
    BigPredictionDf.to_csv('AllSimulations_PredictionSummary.txt',sep='\t')

