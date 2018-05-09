### Transit Tools
### Katelyn Stringer May 6, 2018
### Functions to model the transit of the planet across the star
import numpy as np
#from scipy.optimize import curve_fit


############# FUNCTIONS #######################################

### Read ck solar-type models
def modelReader(temp):
   bstring = np.genfromtxt('./ck_models/ckm00/ckp00_'+str(temp)+'.ascii'\
   ,delimiter='\s+',dtype=None, encoding='utf-8')
   #lstring = [a.decode('utf-8') for a in bstring]
   wave = [float(a.split()[0]) for a in bstring]
   fluxdens = [float(a.split()[-1]) for a in bstring]
   return [0.1*a for a in wave],fluxdens

### Calculate the star's mass from the temperature
def calcStarMass(tstar):
   return 1.989e30*(((tstar/5777)**4)**0.4)

### Calculate the star's radius from the mass
def calcStarRadius(mstar):
   return 6.957e8*(1*(mstar/1.989e30)**0.8)

### Estimate period of orbit
def estimatePeriod(mass_p, mass_s, semi_a):
   return np.sqrt((4*np.pi*np.pi*(semi_a**3))/\
      (6.67e-11*(mass_p + mass_s)))

### Estimate the total transit time
def calcTotalTransitTime(period,rplanet,rstar,semi_a,inclination):
   return (period*rstar)/(np.pi*semi_a)*\
      np.sqrt((1+(rplanet/rstar))**2-\
      ((semi_a/rstar)*np.cos((2*np.pi/360)*inclination))**2)

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
   return (objarea/(np.pi*(rstar**2))) ## cgs units

### Transit wrapper Function
def calcTransit(totflux,period,r_planet,r_star,semia,inclination):
   ### Calculate total transit time
   transtime = calcTotalTransitTime(period,r_planet,r_star,semia,inclination)
   if np.isfinite(transtime)==False:
      transtime = 0

   omega = np.pi/period
   trans_angle = omega*transtime


   times = np.linspace(-2*trans_angle, 2*trans_angle,10000)
   flux = np.ones(len(times))

   ### Convert inclination angle into radians
   inclination = (2*np.pi/360)*inclination
   ### Scale the semimajor axis to the planet's radius
   #semia = semia/r_star
   ### Scale the planet's radius to the planet's radius
   #r_planet = r_planet/r_star
   ### Scale the planet's radius to 1
   #r_star = r_star/r_star

   ### Calculate the impact parameter
   b = semia*np.cos(inclination)/r_star

   ### Determine if the planet is blocking the star
   for j in range(len(times)):
      ### Determine x positions of the planet
      dist = calcDist(semia,times[j],inclination,r_star)
      ### If the planet is not blocking the star
      if (abs(dist) > (r_star+r_planet)):
         pass
      ### If the planet is well past the edge of the star's disk
      elif (abs(dist) <= (r_star - r_planet)):
         area = np.pi*r_planet**2
         flux[j] = 1-calcFluxBlocked(area,r_star)
      ### If the planet is on the edge of the star disk
      else:
         area = calcAreaBlocked(r_planet,r_star,dist)
         flux[j] = 1-calcFluxBlocked(area,r_star)
   return times, times*(period)/(2*np.pi), totflux*flux, trans_angle

### Convolve the flux with the TESS bandpass (only do this once)
def interpTESSFilter(ckwave):
   test = np.genfromtxt('./tess-response-function-v1.0.csv',\
      delimiter=',',skip_header = 8,dtype=None)
   wave = [float(a) for a in test[:,0]]
   throughput = [float(a) for a in test[:,1]]
   tck =  np.interp(ckwave,wave,throughput)
   np.savetxt('./tess_ck_response.csv', np.array([ckwave,tck]).T, \
      delimiter=',')
   return

### Read in the Tess response functions
def read_TESS_response():
   tesstemp = np.loadtxt('./tess_ck_response.csv',delimiter=',')
   return tesstemp[:,0],tesstemp[:,1]

### Multiply the ckflux by the TESS tess_response
def convolveFilter(tess_efficiency,ck_flux,rstar):
   return [a*b for a,b in zip(tess_efficiency,ck_flux)]

### Integrate over flux curve
def integrateFlux(fwave,fcurve):
   return np.trapz(fcurve,x=fwave)

### Calculate magnitudes
def calcMag(fluxes,distance,rstar):
   fluxes = (10**6)*(((rstar)**2)/((3.0857e16*10)**2))*fluxes#distance ### convert pc into meters
   return [-2.5*np.log10(a) - 0.617918 + 5*np.log10(distance)-5 for a in fluxes]#+ 2.5*np.log10(distance**2) for a in fluxes]

### Calculate Vega magnitude
def calcVega():
   vw,vf = modelReader(9500)
   tw_, tr_, tm_, te2_ = setupTransit()
   mv = calcStarMass(9500)
   rv = calcStarRadius(mv)
   vf_ = convolveFilter(tr_,vw,rv)
   tf = integrateFlux(tw_,vf_)
   vegamag = -2.5*np.log10((10**10)*((rv)**2/(3.0857e16*10)**2)*tf)
   return vegamag

### Read in TESS noise function
def readTessError():
   test= np.loadtxt('./tess_noise.txt',delimiter=' ',skiprows=1)
   return test[:,0], test[:,1], test[:,2], test[:,3]
'''
def exponentialFit(A,B,C,x):
   return A*np.exp(B*x)+C

def quadFit(A,B,C,D,E,F,x):
   return A*x**5 + B*x**4 + C*x**3 + D*x**2 + E*x + F

### Fit a function to TESS error curve
def fitError(errmag,noise):
   popt,pcov= curve_fit(quadFit,errmag,np.log10(0.000001*noise))#np.polyfit(errmag, 0.000001*noise, 4)
   return popt
'''
### Determine noise for this magnitude
def calcNoise(errmag,noise,magnitudes):
   return [-2.5*np.log10(a) for a in 0.000001*np.interp(magnitudes,errmag,noise)]#[quadFit(popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],a) \
      #for a in magnitudes]  #popt[0]*a**4 + popt[1]*a**3 + popt[2]*a**2 + popt[3]*a + popt[4]

def calcSNR(magnitudes,magnituderrs):
   return abs(10**(-0.4*(max(magnitudes)-min(magnitudes))))/(10**-0.4*np.mean(magnituderrs))#np.std(magnituderrs)


### Overarching FUNCTIONS
def setupTransit():
   tw, tr = read_TESS_response()
   tm, te1, te2, te27 = readTessError()
   return tw,tr,tm,te2

def doTransit(teff,rp,mp,a,d,i,tw,tr,tm,te2):
   ckw,ckfd = modelReader(teff)
   ms = calcStarMass(teff)
   rs = calcStarRadius(ms)
   #print(ms)
   stflux = convolveFilter(tr,ckfd,rs)
   print(rs)
   ### Integrate over total flux
   tf = integrateFlux(tw,stflux)
   ### Estimate the orbital period
   p = estimatePeriod(ms,mp,a)
   ### Model the transit curve
   phase, t, f, tangle = calcTransit(tf,p,rp,rs,a,i)

   ### Calculate the magnitude
   mags = calcMag(f,d,rs)
   ### Read in TESS errors and fit a function to them
   tm, te1, te2, te27 = readTessError()
   #polys = fitError(tm,te2)
   me = calcNoise(tm,0.000001*te2,mags)
   snr = calcSNR(mags,me)
   ### Calculate the approximate total transit time
   tt = calcTotalTransitTime(p,rp,rs,a,i)
   #print(snr)
   #print(tt)
   #print(tangle)
   #print(phase)
   #print(t)
   #print(mags)
   #print(me)
   #print(p)
   throwaway = (snr,tt,tangle,phase,t,mags,me,p)
   try:
      return snr,tt,tangle,phase,t,mags,me,p
   except:
      return 0,0,0,0,0,0,0,p

################ USAGE ######################################

if __name__== "__main__":

   import matplotlib.pyplot as plt

   ### Define inputs
   ### The example is Sun / Jupiter
   teff = 5750 ### solar effective temperature
   rp = 69.911e6 ### Radius of Jupiter in meters
   mp = 1.898e27 ### mass of Jupiter in kg
   a = (1.496e11)*5.2 ### Semimajor axis of Jupiter in meters
   i = 90#89.9488
   d = 4.8481705933824E-6#10#250#1./(2.06265e5)#10 # distance in parsecs

   ckw,ckfd = modelReader(teff)
   interpTESSFilter(ckw)

   tw, tr = read_TESS_response()

   ### Estimate the stellar mass & radius
   ms = calcStarMass(teff)
   rs = calcStarRadius(ms)

   ### Convolve the filter with the model flux
   stflux = convolveFilter(tr,ckfd,rs)

   ### Integrate over total flux
   tf = integrateFlux(tw,stflux)

   ### Estimate the orbital period
   p = estimatePeriod(ms,mp,a)

   ### Model the transit curve
   phase, t, f, tangle = calcTransit(tf,p,rp,rs,a,i)
   #print(f)
   ### Calculate the magnitude
   mags = calcMag(f,d,rs)

   ### Read in TESS errors and fit a function to them
   tm, te1, te2, te27 = readTessError()
   #polys = fitError(tm,te2)
   me = calcNoise(tm,te2,mags)

   snr = calcSNR(mags,me)


   sn,tt_,ta,ph,t_,ma,merr,per = doTransit(teff,rp,\
               mp,a,d,i,tw,tr,tm,te2)

   ### Calculate the approximate total transit time
   tt = calcTotalTransitTime(p,rp,rs,a,i)
   print('Total transit time = ',str(tt/86400)[0:6],'days')
   print('SNR = ',str(snr)[0:6])

   plt.errorbar([(1./86400)*a - (1./86400)*min(t) for a in t],mags,yerr=me)
   plt.plot([(1./86400)*a - (1./86400)*min(t) for a in t],mags)
   plt.xlabel('Time (days)')
   plt.ylabel('Magnitude (mag)')
   plt.gca().invert_yaxis()
   plt.figtext(0.7,0.5,'Total transit time \n'+str(tt/86400)[0:6]+' days')
   plt.figtext(0.15,0.5,'SNR = '+str(snr)[0:6]+'')
   plt.savefig('./jupiter_transit.png')
   plt.close()


   plt.plot(tw,100*tr)
   plt.tight_layout()
   plt.xlabel('Wavelength (nm)')
   plt.ylabel('Filter Response (%)')
   plt.xlim([400,1250])
   plt.savefig('./tess_filter.png')
   plt.close()

   plt.plot(ckw,ckfd)
   plt.xlabel('Wavelength (nm)')
   plt.tight_layout()
   plt.ylabel('Flux Surface Density (erg/(cm^2 A s))')
   plt.xlim([100,2500])
   plt.savefig('./solar_ck_model.png')
   plt.close()

   plt.figure()
   plt.plot(ckw,stflux)
   plt.tight_layout()
   plt.xlabel('Wavelength (nm)')
   plt.ylabel('Flux Surface Density (erg/(cm^2 A s))')
   plt.xlim([100,1250])
   plt.savefig('./convolved_model.png')
   plt.close()


   #sn,tt_,ta,ph,t_,ma,merr,per = doTransit(14000,637100,mp,1.5e12,166.8101,90,tw,tr,tm,te2)
