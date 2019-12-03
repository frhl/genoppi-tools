# use Anaconda3

# process data from 4 sets os IPs from 4 baits for input into Genoppi

# FOR EACH BAIT:
# (1) log2 transfomation + median normalization of intensity values in each sample
# (2) remove contaminants and proteins with < 2 unique peptides
	# FOR EACH PAIR OF CELL LINES:
	# (3) missing value imputation:
	#     (a) NoImp: remove entries with any missing values across samples
	#	OR
	#     (b) MinImp: replace missing values with randomly sampled values from normal distribution 
	#		  with width 0.3 and down-shifted 1.8 compared to observed sample distribution
	#		  output IMPUTED gene list
	# (4) calcualte logFC: log2(bait) - log2(control) for each bait/control pair
	# (5) generate Genoppi input files 

import pandas as pd
import numpy as np
import sys
sys.path.append('/psych/lagelab/yuhanhsu/Genoppi/src')
from preprocessFunctions import *


# list of baits and corresponding protein reporrts
baitList = ['BCL2','MDM2','PTEN','TARDBP'] # TDP43 = TARDBP (in Uniprot -> gene name mapping talbe)
fileList = ['Genoppi_IPs/'+ x for x in
	['191107LGPSAM06218_RABMONO_BCL2_13-24_HDFNvsAllothers_UniprotHumanNormalizationonFirstSample_Proteins.xlsx',
	'191107LGPSAM06218_RABMONO-MDM2_HDFNvsALL_NormalizedtoFirstSample.xlsx',
	'191107LGPSAM06218_RABMONO-PTEN_HDFNvsALL_NormalizedtoFirstSample.xlsx',
	'191107LGPSAM06218_RABPOLY-Allreruns_HDFNvsALL_NormalizedtoFirstSample.xlsx']]

# list of cell lines
cellList = ['HDFN','G401','T47D','A375']

# iterate through each bait
for i in range(len(baitList)):
	print('##########')
	print(baitList[i])

	# read in raw data
	rawTable = pd.read_excel(fileList[i],sheet_name=0) # read in 1st sheet
	print(rawTable.shape)

	# set "gene" column: map Uniprot ID ('Accession') to gene name (if ID found in ID mapping file)
	rawTable['gene'] = uniprotToGeneName(rawTable['Accession'])
	#rawTable[['gene']].isnull().sum() # make sure no missing value in 'gene' column

	# get relevant columns from MS report and set missing values
	if baitList[i] == 'PTEN':
		intTable = rawTable[['gene','Abundances (Grouped): F15','Abundances (Grouped): F16','Abundances (Grouped): F17',
			'Abundances (Grouped): F18','Abundances (Grouped): F19','Abundances (Grouped): F20',
			'Abundances (Grouped): F21','Abundances (Grouped): F22','Abundances (Grouped): F23',
			'Abundances (Grouped): F24','Abundances (Grouped): F25','Abundances (Grouped): F26']].copy()
	else:
		intTable = rawTable[['gene','Abundances (Grouped): F1','Abundances (Grouped): F2','Abundances (Grouped): F3',
			'Abundances (Grouped): F4','Abundances (Grouped): F5','Abundances (Grouped): F6',
			'Abundances (Grouped): F7','Abundances (Grouped): F8','Abundances (Grouped): F9',
			'Abundances (Grouped): F10','Abundances (Grouped): F11','Abundances (Grouped): F12']].copy()

	intTable.columns = ['gene','HDFN_1','HDFN_2','HDFN_3','G401_1','G401_2','G401_3',
		'T47D_1','T47D_2','T47D_3','A375_1','A375_2','A375_3'] # as indicated in sample submisison form

	#intTable = intTable.replace(0,np.NaN) # replace zeros with NaN
	print(intTable.shape)


	# (1) log2 transfomation + median normalization of intensity values in each sample
	intTable.iloc[:,1:] = np.log2(intTable.iloc[:,1:])
	normTable = medianNorm(intTable)
	print(normTable.shape)

	# (2) remove contaminants and proteins with < 2 unique peptides
	filterInds = (~rawTable['Contaminant']) & (rawTable['# Unique Peptides'] >= 2) # rawTable and normTable entries are in same order
	filterTable = normTable[filterInds]
	print(filterTable.shape)


	# iterate through eac pair of cell lines
	# (3) missing value imputation:
	#     (a) NoImp: remove entries with any missing values across samples
	#	OR
	#     (b) MinImp: replace missing values with randomly sampled values from normal distribution 
	#		  with width 0.5 and down-shifted 1.8 compared to observed sample distribution
	# (4) calcualte logFC: log2(bait) - log2(control) for each bait/control pair
	# (5) generate Genoppi input files 

	for j in range(0,len(cellList)-1):
		for k in range(j+1,len(cellList)):

			print('#---------')
			print(cellList[j] + ' vs. ' + cellList[k])
		
			pairTable = filterTable[['gene',cellList[j]+'_1',cellList[k]+'_1',
				cellList[j]+'_2',cellList[k]+'_2',cellList[j]+'_3',cellList[k]+'_3']]

			# NoImp
			noImpTable = pairTable.dropna() # drop entries with any missing values
			noImpTable = calcLogFC_fromLogValues(noImpTable)
			print(noImpTable.shape)
			noImpTable.to_csv('GenoppiInput/'+baitList[i]+'.'+cellList[j]+'vs'+cellList[k]+
				'.NoImp.GenoppiInput.txt',sep="\t",index=False)

			# MinImp
			minImpTable = pairTable[pd.notnull(pairTable).sum(axis=1) > 0] # drop entries with all missing values
			np.random.seed(123) # so rerunning gives same imputation results
			minImpTable = minNormImpute(minImpTable,1.8,0.3) # down-shifted by 1.8*SD, width of 0.3*SD

			# ouput "IMPUTED" gene list 
			impGeneList = minImpTable.iloc[:,[0,-1]]
			impGeneList.columns = ['gene','significant']
			impGeneList.to_csv('GenoppiInput/'+baitList[i]+'.'+cellList[j]+'vs'+cellList[k]+'.ImpGeneList.txt',
				sep="\t",index=False)

			minImpTable = calcLogFC_fromLogValues(minImpTable.iloc[:,:-1])
			print(minImpTable.shape)
			minImpTable.to_csv('GenoppiInput/'+baitList[i]+'.'+cellList[j]+'vs'+cellList[k]+
				'.MinImp.GenoppiInput.txt',sep="\t",index=False)

print('SCRIPT COMPLETED')