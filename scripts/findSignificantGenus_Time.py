import numpy as np
import pandas as pd
import scipy.stats as stats

P_CUTOFF = .050

## LOAD DATA
ClassResults = pd.read_csv('../data/take-home_classification_results.csv',index_col='classification_id')
MetaData=pd.read_csv('../data/take-home_metadata.csv')
TaxaInfo=pd.read_csv('../data/take-home_tax_info.csv')

## Replace id with genus name
ClassResults.columns=TaxaInfo.name

geni = ClassResults.columns.values.tolist()


## FIND Taconic mice gavaged with JAX and time data

TAC_galJAX=MetaData.loc[MetaData['isolation_source_s'].str.contains("Taconic mice gavaged with JAX ")]

# Isolate time data from source column
timeCapturePat = '(\d+)([dh]$)'
TAC_galJAX[['time','timeUnit']]=TAC_galJAX['isolation_source_s'].str.extract(timeCapturePat)
TAC_galJAX['time'] = pd.to_numeric(TAC_galJAX.time, errors='coerce')

# Replace hours with days to keep units consistent, remove units column (done since only hours and days detected as potential units)
TAC_galJAX['time'] = np.where(TAC_galJAX['timeUnit'] =='h',TAC_galJAX['time']/24,TAC_galJAX['time'])
TAC_galJAX=TAC_galJAX.drop(columns='timeUnit')

TAC_galJAX=TAC_galJAX.set_index('classification_id')


Class_TAC_galJAX = TAC_galJAX.join(ClassResults)


## Finding significant change
sigList=[]
sigVals=[]
for genus in geni:
    relData= Class_TAC_galJAX[['time',genus]]
    timeIncs=list(set(relData['time']))
    testData=[]
    for t in timeIncs:
        tempData= relData.loc[relData['time'] == t]
        testData.append(tempData[genus].values.tolist())

    [f,p]=stats.f_oneway(testData[0],testData[1],testData[2],testData[3],testData[4])
    if p < P_CUTOFF:
        sigList.append(genus)
        sigVals.append(p)

sigDF=pd.DataFrame({'Genus':sigList,'P-Value':sigVals})

# Saving Significant Data
Long_Class_TAC_galJAX = Class_TAC_galJAX.melt(id_vars='time',value_vars=geni,var_name='Genus',value_name='Read Counts')


Sig_Long_Class_TAC_galJAX = Long_Class_TAC_galJAX.merge(sigDF,how='inner', on='Genus')

Sig_Long_Class_TAC_galJAX.to_csv('../results/SignificantGenus_P05_TimeChange_Reads.csv',index=False)
sigDF.sort_values('P-Value').to_csv('../results/SignificantGenus_P05_TimeChange.csv',index=False)

