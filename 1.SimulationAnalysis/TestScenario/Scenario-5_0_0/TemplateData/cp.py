import os,glob

for FileName in glob.glob('../Iter*'):

	print('cp 001.Output_fasta_file_ver4.py {0}'.format(FileName))
	os.system('cp 001.Output_fasta_file_ver4.py {0}'.format(FileName))
