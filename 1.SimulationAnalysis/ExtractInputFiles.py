import os, glob

if __name__ == '__main__':

	MutRateList = sorted(glob.glob('MutRate_*'))
	for MutRate in MutRateList:

		ScenarioList = sorted(glob.glob(MutRate+'/Scenario-*'))
		for Scenario in ScenarioList:

			os.chdir(Scenario)

			Commandline = 'tar -zxvf SimData.tar.gz'
			print('Extracting input files in: %s' % Scenario);os.system(Commandline)

			os.chdir('../../')
