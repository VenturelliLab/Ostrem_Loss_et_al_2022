import numpy as np
import pandas as pd
from scipy.integrate import odeint


def res_cons(y, t, w, C, m, num_res, num_species):
    
    dydt = np.zeros(num_res + num_species)
    # Constructing the equations for species equations
    for sp_in in range(0, num_species):
        for r_in in range(0, num_res):
            dydt[sp_in] = dydt[sp_in] + w[r_in, sp_in]*C[r_in, sp_in]*y[num_species + r_in]
            # i+r_in denotes the resource equation index
        dydt[sp_in] = y[sp_in]*(dydt[sp_in] - m[sp_in])
        
    # Constructing the equations for resource requations.
    for r_in in range (0, num_res):
        for sp_in in range(0, num_species):
            dydt[r_in + num_species] = dydt[r_in + num_species] + -1*C[r_in, sp_in]*y[r_in + num_species]*y[sp_in]
    return dydt



MasterDF=pd.read_csv('MasterDF.csv')
subDF=MasterDF[(MasterDF['Passage']==0)&(MasterDF['Media']=='Inulin')]
comm10species=['AC','DP','BH','DL','EL','BT','BV','PJ','BP','CS']
conlist=['COMM10<PJ','COMM10<BT','COMM10<BV','COMM10<DP','COMM10<BP','COMM10<EL','COMM10<BH','COMM10<AC','COMM10<CS','COMM10<DL','COMM10<equal']
num_res = 1 # Total number of resources
num_species = 10 # Total number of species
 # Resource uptake matrix
m = np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]) # Maintenance vector
  
df=pd.DataFrame()
z=0
for rep in range(1,101):

    for numuts in range(1,11):
        for gparam in range(2,16):
            print(str(gparam*0.1))
            w1=np.random.normal(gparam*0.1, .1, size=10)#change variance, sample variances
            C1=np.random.normal(gparam*0.1, .1, size=10)
            w1[numuts:11]=0
            C1[numuts:11]=0
            w=np.zeros((1,10))
            w[0]=w1
            C=np.zeros((1,10))
            C[0]=C1

            res1=1
            init_vals_res = np.array([res1]) # initial values for resources

            t = np.linspace(0, 24, 50)
            for treatment in conlist: 

                tempDF=subDF[subDF['Treatment']==treatment]
                abund=[]
                for species in comm10species:
                    abund.append(list(tempDF[species+' Absolute Abundance'])[0])
                init_vals_species = np.array(abund)
                init_vals = np.hstack((init_vals_species, init_vals_res))
                sim_data = odeint(res_cons, init_vals, t, args = (w, C, m, num_res, num_species))


                #Storing the data
                df.at[z,'Number of Utilizers']=numuts
                totalOD=[]
                for sp in range(num_species):
                    df.at[z,'Species '+str(sp)+' Abundance']=sim_data[49,sp]
                    totalOD.append(sim_data[49,sp])
                totalOD=sum(totalOD)
                for sp in range(num_species):
                    df.at[z,'Species '+str(sp)+' Fraction']=sim_data[49,sp]/totalOD
                df.at[z,'Species with high initial abundance']=treatment
                df.at[z,'Rep']=rep
                df.at[z,'Resource']=res1
                df.at[z,'Median Growth']=gparam*0.1

                z+=1
df.to_csv('Simulations.csv')