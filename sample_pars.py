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

### Get 9 planet sizes
r_p = 6.371e6*np.logspace(-1,3,1)#9)

### Get 5 semimajor axis values
au = (1.496e11)*np.logspace(-1,2,5)

### Get 10 distances
d_s = np.logspace(0,2.5,1)#0)

### Get total number of iterations
tnum = len(temps)*len(r_p)*len(au)*len(d_s)

### Set up pandas dataframe for all outcomes
stardata = pd.DataFrame({'temp':np.zeros(tnum),'rplanet':np.zeros(tnum),\
   'dist':np.zeros(tnum),'a':np.zeros(tnum),'snr':np.zeros(tnum),'transit_time':np.zeros(tnum),\
   'period':np.zeros(tnum),'transitangle':np.zeros(tnum)})

mpl = 1.898e27
inc = 90



### Calculate the transit curve for each grid point
for t in range(len(temps)):
   for r in range(len(r_p)):
      for d_ in range(len(d_s)):
         for a_ in range(len(au)):
            stardata['temp'].iloc[t+r+d_+a_] = temps[t]
            stardata['rplanet'].iloc[t+r+d_+a_] = r_p[r]
            stardata['dist'].iloc[t+r+d_+a_] = d_s[d_]
            stardata['a'].iloc[t+r+d_+a_] = au[a_]

            sn,tt_,ta,ph,t_,ma,merr,per = doTransit(temps[t],r_p[r],\
               mpl,au,d_s[d_],inc,tw_,tr_,tm_,te2_)

            stardata['snr'].iloc[t+r+d_+a_] = sn
            stardata['transit_time'].iloc[t+r+d_+a_] = tt_
            stardata['period'].iloc[t+r+d_+a_] = per
            stardata['transitangle'].iloc[t+r+d_+a_] = ta

            lc = pd.DataFrame({'time':ph,'phase':t_,'mag':ma,'magerr':merr})
            lc.to_csv('./output_lc/'+str(int(temps[t]))+'_'+str(r_p[r])+\
               '_'+str(d_s[d_])+'_'+str(au[a_])+'.csv',header=True,index=False)

stardata.to_csv('./star_results.csv',header=True,index=False)
