import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from transit_tools import *

# Plot different slices in parameter space

### Read in files
sd = pd.read_csv('./star_results.csv',header=0)

### Remove crazy valyes
sd = sd.loc[sd['snr']>0]#sd.loc[((sd['snr'])<200)&(sd['snr']>0)]
#sd = sd.loc[(sd['rplanet'])<200]

ssd = sd.loc[sd['dist']==1]

sds = sd.loc[sd['temp']==5750]

### Plot rplanet vs temp of star, w/ snr colors
plt.figure(figsize=(12,10))
plt.scatter(ssd['temp'],np.log10((1./6.371e6)*(ssd['rplanet'])),\
   c=[1/a for a in ssd['snr']])
plt.xlabel('Stellar Effective Temperature (K)')
plt.ylabel('Exoplanet Radius (R_Earth)')
#plt.ylim([0,20])
plt.colorbar()
plt.figtext(0.85,0.5,'Signal-to-Noise Ratio',rotation=270)
plt.title('Fixed Distance = 1 pc')
plt.savefig('./t_r_snr.png')
plt.close()


plt.figure(figsize=(12,10))
plt.scatter(sds['dist'],np.log10((1./6.371e6)*(sds['rplanet'])),\
   c=[(1/a) for a in sds['snr']])
plt.xlabel('Distance (pc)')
plt.ylabel('Exoplanet Radius (R_Earth)')
#plt.ylim([0,20])
plt.colorbar()
plt.figtext(0.85,0.5,'Signal-to-Noise Ratio',rotation=270)
plt.title('Fixed Temperature = 5750')
plt.savefig('./d_r_snr.png')
plt.close()

plt.hist(sd['dist'],bins=100)
plt.show()
