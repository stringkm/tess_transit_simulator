import numpy as np

def read_ck_ascii(filename):
   lamb,fluxdens = np.loadtxt(filename,sep='\s+',header=None)
   return lamb,fluxdens
