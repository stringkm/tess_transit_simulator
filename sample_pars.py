import numpy as np
from transit_tools import *
from glob import glob
import pandas as pd
import time

### Import the basics to set up the transit code
tw_, tr_, tm_, te2_ = setupTransit()

st = time.time()

##### Create grid of parameters

### Read in the solar metallicity star files
fnames = glob('./ck_models/ckm00/*.ascii')

### Get a list of all stellar types in folder
temps = [int(str(a).split('_')[-1].split('.')[0]) for a in fnames]

### Get 10 planet sizes
r_p = [6.371e6*10**a for a in np.linspace(-1,2.5,10)]#0)]#9)]#9)

### Get 5 semimajor axis values
au = [(1.496e11)*10**a for a in [1]]#np.logspace(-1,1,5)]

### Get 10 distances
d_s = [10**a for a in np.linspace(0,2.5,10)]#10)]#0)

### Get total number of iterations
tnum = len(temps)*len(r_p)*len(au)*len(d_s)
print(tnum)

### Set up pandas dataframe for all outcomes
stardata = pd.DataFrame({'temp':np.zeros(tnum),'rplanet':np.zeros(tnum),\
   'dist':np.zeros(tnum),'a':np.zeros(tnum),'snr':np.zeros(tnum),'transit_time':np.zeros(tnum),\
   'period':np.zeros(tnum),'transitangle':np.zeros(tnum)})

mpl = 1.898e27
inc = 90

iterator = 0

### Calculate the transit curve for each grid point
for t in range(len(temps)):
   for r in range(len(r_p)):
      for d_ in range(len(d_s)):
         for a_ in range(len(au)):

            stardata['temp'].iloc[iterator] = temps[t]
            stardata['rplanet'].iloc[iterator] = r_p[r]
            stardata['dist'].iloc[iterator] = d_s[d_]
            stardata['a'].iloc[iterator] = au[a_]

            if 1==1:
               sn,tt_,ta,ph,t_,ma,merr,per = \
                  doTransit(temps[t],r_p[r],mpl,au[a_],d_s[d_],inc,tw_,tr_,tm_,te2_)

               stardata['snr'].iloc[iterator] = sn
               stardata['transit_time'].iloc[iterator] = tt_
               stardata['period'].iloc[iterator] = per
               stardata['transitangle'].iloc[iterator] = ta

               lc = pd.DataFrame({'time':ph,'phase':t_,'mag':ma,'magerr':merr})
               lc.to_csv('./output_lc/'+str(int(temps[t]))+'_'+str(r_p[r])+\
                   '_'+str(d_s[d_])+'_'+str(au[a_])+'.csv',header=True,index=False)
               iterator = iterator + 1
               print(time.time()-st,'s elapsed')
               st = time.time()
            else:
               print(iterator,'failed')
               iterator = iterator + 1
               pass
stardata.to_csv('./star_results.csv',header=True,index=False)
