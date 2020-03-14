import os, glob

Release_version = '93'

fin = open('download.list','r')
for line in fin:
	species = line.strip()

	## genome ##
	code = 'wget ftp://ftp.ensembl.org/pub/release-'+Release_version+'/fasta/'+species+'/dna/*.dna.toplevel.fa.gz'

	print(code);os.system(code)


	## cds ##
	code = 'wget ftp://ftp.ensembl.org/pub/release-'+Release_version+'/fasta/'+species+'/cds/*.cds.all.fa.gz'

	#print(code);os.system(code)


	## pep ##
	code = 'wget ftp://ftp.ensembl.org/pub/release-'+Release_version+'/fasta/'+species+'/pep/*.pep.all.fa.gz'

	#print(code);os.system(code)


	## gtf ##
	code = 'wget ftp://ftp.ensembl.org/pub/release-'+Release_version+'/gtf/'+species+'/*.'+Release_version+'.gtf.gz'

	print(code);os.system(code)

	
os.system('gunzip *.gz')
os.system('ls | while read upName; do loName=`echo "${upName}" | tr \'[:upper:]\' \'[:lower:]\'`; mv "$upName" "$loName"; done')
