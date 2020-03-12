import os,glob

## This code is used to remove unnecessary files generated in the simulation

fin = open('RemoveFiles.list','r');FileList = []
for line in fin:FileList.append(line.strip())
fin.close()

Dir1List = sorted(glob.glob('Scenario-*'))
for Dir1 in Dir1List:

	print('Removing unnecessary files in: %s' % Dir1)
	os.chdir(Dir1)
	for FileName in FileList:
		Code = "find -maxdepth 5 -name '%s' -exec rm {} \;" % FileName
		os.system(Code)

	Dir2List = sorted(glob.glob('Iter*'))
	for Dir2 in Dir2List:
		
		VcfList = sorted(glob.glob(Dir2+'/*.vcf'))
		for Vcf in VcfList:
			if Vcf.replace('.vcf','').split('/')[-1] in os.listdir(Dir2):
				os.system('rm -r '+Vcf.replace('.vcf',''))

	os.chdir('../')
