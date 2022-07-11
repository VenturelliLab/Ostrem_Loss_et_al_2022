import pandas as pd
import numpy as np
from scipy.spatial import distance
predDF=pd.read_csv('Simulations_AllNegative.csv')
distDF0=pd.DataFrame()
species = ['S'+str(i+1) for i in range(10)]
num2species={
        8: 'CS',
5:'DL',
6: 'BH',
2: 'AC',
7: 'BP',
1: 'BV',
0: 'PJ',
4:'BT',
9: 'EL',
3: 'DP'}
k=0

for rep in range(100):
    
    for med in [-1,-0.5,-0.25,-0.1,0]:

        for passage in range(4):
            print(passage)
            subdf=predDF[(predDF['Median']==med)&(predDF['Passage']==passage)&(predDF['Rep']==rep)]

            blahDF=subdf.dropna()
            newdf=blahDF.reset_index(drop=True)
            newDF=newdf.drop(columns=['Passage','Treatment','Median','Rep','% Positive','% Negative'])
            for spe in range(len(species)):
                newDF=newDF.drop(columns=[num2species[spe]+' Absolute Abundance'])
            trep=0
            print(newDF)
            for z in list(newDF.index.values):
                for i in list(newDF.index.values):
                    if i>z:
                        x=np.array(newDF.loc[z])
                        y=np.array(newDF.loc[i])
                        dist=distance.euclidean(x,y)
                        distDF0.at[k,'Euclidean Distance']=dist
                        distDF0.at[k,'Bio Rep']=rep
                        distDF0.at[k,'Tech Rep']=trep
                        distDF0.at[k,'Median']=med
                        distDF0.at[k,'Passage']=passage
                        k+=1 
                        trep+=1
                    else:
                        continue
distDF0.to_csv('EuclideanDistances_AllNegative_moremeds.csv')
