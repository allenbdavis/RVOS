import numpy as np
from RVOS_Constants import *

class Star:
    def __init__(self,massStar,t0,coord,sigma,planets=None,name=''):
        self.typecheck = 'Star'
        self.massStar = massStar #in mSun. Default: 1 solar mass
        self.name = name
        self.t0 = t0 #nominal start date of observations
        self.coord = coord #ra/dec position on the sky
        self.sigma = sigma #uncertainty in RV measurements
        if not planets == None:
            self.planets = planets
        
class Planet:
    
    def __init__(self,star,mass,period,ecc,incl,w,tp,K=None,isReal=True):
        self.typecheck = 'Planet'
        self.isReal = isReal
        self.star = star
        self.name = star.name #initialized to star's name. Will append letter later.
        self.mass = mass #in mEarth
        self.period = period #in years
        self.massStar = star.massStar #in mSun. Default: 1 solar mass
        self.ecc = ecc #eccentricity. Default: 0.2
        self.incl = incl #inclination in degrees. Default: 90.
        self.sini = np.sin(incl*np.pi/180.) #sin(i)
        self.t0 = star.t0 #JD of start. default: "today"

        def getSMA():
            totmass = (self.mass*mEarth)+(star.massStar*mSun)
            return (1/au) * ( ( (self.period*year)**2 * G*totmass ) / ( 4.*pi*pi))**(1./3.)
        
        self.sma = getSMA() #semimajor axis in AU. Given uniquely by mass, massStar, and period.
        
        if w==None:
            self.w = 2.*pi*np.random.random() #phase. default: random phase from [0,2pi)
        else:
            self.w = w
        
        if tp==None:
            self.tp = self.t0+(np.random.random()*self.period) #periastron time. Default: random from [today,today+period)
        else:
            self.tp = tp
            
        if K == None:
            self.K = (2.*pi*G/(self.period*day))**(1./3.) * (self.mass*mEarth*self.sini) * ((self.massStar*mSun)**(-2./3.)) * (1./np.sqrt(1.-self.ecc**2.))
        else:
            self.K = K #used for a fitted planet, where the K is known but the mass is actually just msini
        self.params = [self.period,self.ecc,self.tp,self.w,self.K,self.mass]
        self.params_print = [self.period,self.ecc,self.tp-self.t0,self.w*180./pi,self.K,self.mass]
        
