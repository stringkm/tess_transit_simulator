### Transit Tools
### Katelyn Stringer May 6, 2018
### Functions to model the transit of the planet across the star
import numpy as np

############# FUNCTIONS #######################################

### Calculate the star's mass from the temperature
def calcStarMass(tstar):
   return 1.989e30*(((tstar/5777)**4)**0.4)

### Calculate the star's radius from the mass
def calcStarRadius(mstar):
   return 6.957e8*(1*(mstar/1.989e30)**0.8)

### Estimate period of orbit
def estimatePeriod(mass_p, mass_s, semi_a):
   return np.sqrt((4*np.pi*np.pi*(semi_a**3))/\
      (6.67e-11*(mass_p+mass_s)))

### Calculate position for each update
def calcDist(semi_a,objtime,obji,rstar):
   dd = semi_a*np.sqrt((np.sin(objtime)**2)+\
      (np.cos(obji)**2)*(np.cos(objtime)**2))
   return dd

### Calculate the area of star blocked by the transiting planet
def calcAreaBlocked(rplanet,rstar,dd):
   p1 = (rplanet**2)*np.arccos((dd**2 + rplanet**2 - \
      rstar**2)/(2*rplanet*dd))
   p2 = (rstar**2)*np.arccos((dd**2 + rstar**2 - \
      rplanet**2)/(2*rstar*dd))
   p3 = -0.5*np.sqrt((-1*dd+rplanet+rstar)*(dd+rplanet-rstar)\
      *(dd-rplanet+rstar)*(dd+rplanet+rstar))
   return p1+p2+p3

### Calculate the amount of star flux blocked by planet
def calcFluxBlocked(objarea,rstar):
   return (objarea/(np.pi*rstar**2))

### Transit wrapper Function
def calcTransit(r_planet,r_star,semia,inclination):
   times = np.linspace(-0.05*np.pi, 0.05*np.pi,100)
   flux = np.ones(len(times))
   ### Convert inclination angle into radians
   inclination = (2*np.pi/360)*inclination
   ### Scale the semimajor axis to the planet's radius
   semia = semia/r_star
   ### Scale the planet's radius to the planet's radius
   r_planet = r_planet/r_star
   ### Scale the planet's radius to 1
   r_star = r_star/r_star

   ### Calculate the impact parameter
   b = semia*np.cos(inclination)/r_star

   print(r_star,r_planet,semia)
   ### Determine if the planet is blocking the star
   for j in range(len(times)):
      ### Determine x positions of the planet
      dist = calcDist(semia,times[j],inclination,r_star)
      ### If the planet is not blocking the star
      if (abs(dist) > r_star):
         pass
      ### If the planet is well past the edge of the star's disk
      elif (abs(dist) < (r_star-r_planet)):
         area = np.pi*r_planet**2
         flux[j] = 1-calcFluxBlocked(area,r_star)
      ### If the planet is on the edge of the star disk
      else:
         print(dist/r_planet)
         area = calcAreaBlocked(r_planet,r_star,dist)
         flux[j] = 1-calcFluxBlocked(area,r_star)
   return times, flux

################ USAGE ######################################

if __name__== "__main__":
   ### Define inputs
   ### The example is Sun / Jupiter
   teff_s = 5777 ### solar effective temperature
   rp = 69.911e6 ### Radius of Jupiter in meters
   mp = 1.898e27 ### mass of Jupiter in kg
   a = (1.496e11)*0.1 ### Semimajor axis of Jupiter in meters

   ### Estimate the stellar mass & radius
   ms = calcStarMass(teff_s)
   rs = calcStarRadius(ms)

   ### Estimate the orbital period
   p = estimatePeriod(ms,mp,a)

   ### Model the transit curve
   t,f = calcTransit(rp,rs,a,90)
