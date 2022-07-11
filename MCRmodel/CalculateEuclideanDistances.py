import scipy as sp
from scipy.spatial import distance
import pandas as pd
import numpy as np
num_species=10
df=pd.read_csv('Simulations.csv')
distDF=pd.DataFrame()
k=0
for rep in range(1,101):
    for numuts in range(1,11):
        for gparam in [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]:
                        
            subdf=df[(df['Number of Utilizers']==numuts)&(df['Rep']==rep)&(df['Median Growth']==gparam)]

            blahDF=subdf.dropna()
            newdf=blahDF.reset_index(drop=True)
            newDF=newdf.drop(columns=['Number of Utilizers','Species with high initial abundance','Rep','Resource','Median Growth'])
            
            for sp in range(num_species):
                newDF=newDF.drop(columns=['Species '+str(sp)+' Abundance'])
            print('median growth= '+str(gparam))
            print(newDF)

            for z in list(newDF.index.values):
                for i in list(newDF.index.values):
                    if i>z:

                            x=np.array(newDF.loc[z])
                            y=np.array(newDF.loc[i])
                            dist=distance.euclidean(x,y)
                            distDF.at[k,'Distance Measurent']='Relative'
                            distDF.at[k,'Euclidean Distance']=dist
                            distDF.at[k,'Number of Utilizers2']='m'+str(numuts)
                            distDF.at[k,'Number of Utilizers']=numuts
                            distDF.at[k,'Rep']=rep
                            distDF.at[k,'Median Growth']=gparam
                            k+=1 
                    else:
                        continue
distDF.to_csv('EuclideanDistances.csv')