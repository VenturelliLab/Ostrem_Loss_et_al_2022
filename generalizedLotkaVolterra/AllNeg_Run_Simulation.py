import pandas as pd
import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import minimize, minimize_scalar

def system(x, r, Ax):
    # dx / dt = x (r + A @ x)
    return x*(r + Ax)

def runODE(t_eval, x, params):
    # reshape parameters to matrices
    dimx  = len(x)
    r  = params[:dimx]
    A  = np.reshape(params[dimx:], [dimx, dimx])

    # x will be 1-D
    def dX_dt(t, x):
        # compute matrix vector products
        Ax = A@x
        return system(x, r, Ax)

    # solve ODE model
    soln = solve_ivp(dX_dt, (0, t_eval[-1]), x,
                     t_eval=t_eval, method='LSODA')
    t_solver, y = soln.t, soln.y.T
    return y

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

species = ['S'+str(i+1) for i in range(10)]

iDF=pd.read_csv('MasterDF.csv')
conlist=[
'COMM10<AC',
 'COMM10<BH',
 'COMM10<DP',
 'COMM10<BT',
 'COMM10<DL',
 'COMM10<EL',
 'COMM10<equal',
 'COMM10<BP',
 'COMM10<BV',
 'COMM10<PJ']
predDF=pd.DataFrame(columns=['Median','Passage','Treatment','Rep'])
z=0
for rep in range(100):
    mu=np.random.normal(0.25,0.2,size=110)
    aii=np.random.normal(-1.25,0.5,size=110)
    
    for med in [-1,-0.5,-0.25,-0.1,0]:

        simparams=np.random.normal(med, .25, size=110)
        for i in range(10):
            simparams[i]=mu[i]
        for i in range(10):
            simparams[i*11+10]=aii[i*11+10]
        
        for i in range(len(simparams)):
            if simparams[i]>0:
                simparams[i]=0

        estsimparams=list(simparams)
        del estsimparams[0:10]
        del estsimparams[0::11]
        perpos= len([i for i in estsimparams if i >= 0])
        perneg= len([i for i in estsimparams if i < 0])

        for treatment in conlist:

 
            for passage in range(4):
                print(passage)
                if passage==0:
                    subDF=iDF[(iDF['Media']=='Inulin')&(iDF['Passage']==0)&(iDF['Treatment']==treatment)]
                    x0 = np.zeros(len(species))
                    for spe in range(len(species)):
                        x0[spe]=(subDF[num2species[spe]+' Absolute Abundance'])
                else:
                    subDF=predDF[(predDF['Median']==med)&(predDF['Passage']==passage-1)&(predDF['Treatment']==treatment)&(predDF['Rep']==rep)]
                    x0 = np.zeros(len(species))
                    for spe in range(len(species)):
                        x0[spe]=(subDF[num2species[spe]+' Relative Abundance'])/20


                t_eval = np.linspace(0, 24, 50)
                y = runODE(t_eval, x0, simparams)

                ods=[]
                for spe in range(len(species)):
                    predDF.at[z,num2species[spe]+' Absolute Abundance']=y[-1][spe]
                    ods.append(y[-1][spe])
                total=sum(ods)
                for spe in range(len(species)):
                    predDF.at[z,num2species[spe]+' Relative Abundance']=y[-1][spe]/total
                predDF.at[z,'Treatment']=treatment
                predDF.at[z,'Median']=med
                predDF.at[z,'Passage']=passage
                predDF.at[z,'Rep']=rep
                predDF.at[z,'% Positive']=perpos/(perpos+perneg)
                predDF.at[z,'% Negative']=perneg/(perpos+perneg)
                z+=1
                predDF.to_csv('Simulations_AllNegative.csv')
