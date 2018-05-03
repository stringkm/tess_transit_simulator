### Transit function
import numpy as np
from scipy.intergrate import odeint
import argparse

if __name__== "__main__":
   rs = 1
   rp = (7.149e7)/rs
   d = (1.496e11*(1./rs))*[0.1,0.5,1.]

#else:
#  parser = argparse.ArgumentParser()

## Initialize time array
### Estimate period of orbit

def estimatePeriod(mass_p, mass_s, semi_a):
   return np.sqrt((4*np.pi*np.pi*(semi_a**3))/\
      (6.67e-11*(mass_p+mass_s)))


### Initialize time array
t = np.linspace(-0.05*np.pi, 0.05*np.pi,1000)

### Calculate position for each ode update
def calcX(semi_a,objtime,obji,rstar):
   dd = semi_a*np.sqrt((np.sin(objtime)**2)+\
      (np.cos(obji)**2)*(np.cos(objtime)**2))
   return dd-rstar

### Calculate the area of star blocked by the transiting planet
def calcArea(rplanet,xx):
   return (rplanet**2)*np.arccos(xx/rplanet)-\
      rplanet*xx*np.sqrt(1.-(xx/rplanet)**2)

### Calculate the amount of star flux blocked by planet
def calcFluxBlocked(objarea,rstar):
   return (objarea/(np.pi*rstar**2))

### Subtract blocked flux from the star's BB curve
