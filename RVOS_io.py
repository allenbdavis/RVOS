import sys
from RVOS_Bodies import *
from astropy.coordinates import SkyCoord
from astropy import units as u

def saveSystem(star,planets,fname=''):
	from time import strftime
	import os
	import errno
	s=''
	s+='#This file contains the info necessary to recreate the 	Star, Planets, and observation data.\n'
	s+='\nSTAR\n'
	s+='name=%s\n'%star.name
	s+='massStar=%s\n'%star.massStar
	s+='t0=%s\n'%star.t0
	s+='ra=%s\n'%star.coord.ra
	s+='dec=%s\n'%star.coord.dec
	s+='nPlanets=%s\n'%star.nPlanets
	s+='nPlanets_fit=%s\n'%star.nPlanets_fit
	s+='sigma=%s\n'%star.sigma
	
	for i in range(0,star.nPlanets):
		pl = star.planets[i]
		s+='\nREAL PLANET %d\n'%(i+1)
		s+='name=%s\n'%pl.name
		s+='mass=%s\n'%pl.mass
		s+='period=%s\n'%pl.period
		s+='ecc=%s\n'%pl.ecc
		s+='incl=%s\n'%pl.incl
		s+='w=%s\n'%pl.w
		s+='tp=%s\n'%pl.tp
		s+='k=%s\n'%pl.K
	
	for i in range(0,star.nPlanets_fit):
		pl = star.planets_fit[i]
		s+='\nFIT PLANET %d\n'%(i+1)
		s+='name=%s\n'%pl.name
		s+='mass=%s\n'%pl.mass
		s+='period=%s\n'%pl.period
		s+='ecc=%s\n'%pl.ecc
		s+='incl=%s\n'%pl.incl
		s+='w=%s\n'%pl.w
		s+='tp=%s\n'%pl.tp
		s+='k=%s\n'%pl.K
		
		try:
			dum = pl.uncs
			s+='UNCS EXIST\n'
			s+='uncs[5]=%s\n'%pl.uncs[5] #mass
			s+='uncs[0]=%s\n'%pl.uncs[0] #period
			s+='uncs[1]=%s\n'%pl.uncs[1] #ecc
			s+='uncs[3]=%s\n'%pl.uncs[3] #w
			s+='uncs[2]=%s\n'%pl.uncs[2] #tp
			s+='uncs[4]=%s\n'%pl.uncs[4] #k
		except AttributeError:
			pass
	
	s+='\nOBSERVATIONS\n'
	s+='#JD,RV_True,RV_Obs\n'
	for i in range(0,len(star.obsList)):
		s+='%s,%s,%s\n'%(star.obsList[i],star.RV_True[i],star.RV_Obs[i])

	# Make a new directory every day
	path = '.\\figures_'+strftime("%Y-%m-%d")
	def make_sure_path_exists(path):
		try:
			os.makedirs(path)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise
	make_sure_path_exists(path)
	if not fname=='': fname+='_'
	f = open(fname+strftime("%H-%M-%S")+'.txt','w')
	f.write(s)
	f.close()

def loadSystem(fname):
	f = open(fname,'r')
	s = f.readlines()
	f.close()
	
	def isParameter(line):
	#lines start with a lower case letter IFF they are parameter lines
		if line=='':
			return False
		else:
			return line[0].islower()

	def isNewObject(line):
		#lines start with an upper case letter IFF they are new object lienes
		if line=='':
			return False
		else:
			return line[0].isupper()
		
	def isObservation(line):
		#lines start with a numeral IFF they are observation data
		if line=='':
			return False
		else:
			return line[0].isdigit()

	def isRaDec(line):
		#RA and DEC must be handled differently to convert to pure degrees
		return line.startswith('ra=') or line.startswith('dec=')
		
	def isName(line):
		#names must be stored as strings
		return line.startswith('name')
	
	def tryToStore(arr,ind,val):
	# Recursively attempts to store val in arr at ind.
	# This method means that Planets can be specified
	# in any order and without knowing nPlanets a priori.
		try:
			arr[ind] = val
		except IndexError:
			arr.append(None)
			tryToStore(arr,ind,val)
	
	# Will hold all the parameters needed to make a fully realized star
	class ProtoStar:
		def __init__(self):
			pass
	# Will hold all the parameters needed to make a fully realized planet
	class ProtoPlanet:
		def __init__(self):
			pass
	
	# Initialize arrays to hold parameters
	obj = None #current object. Will either be ProtoStar or ProtoPlanet
	ps = ProtoStar()
	pp_arr = [] #array to hold ProtoPlanets()
	pp_fit_arr = [] #array to hold ProtoPlanets() for fitted planets
	obsList = []
	RV_True = []
	RV_Obs = []
	
	# Read file
	for line_raw in s:
		line = line_raw.rstrip()
		# execute a parameter assignment
		if isParameter(line):
			assert not obj==None,'OBJECT must be specified before its parameters.'
			spl1 = line.split('#') #split out possible comments
			spl2 = spl1[0].split('=') #split apart assignment statement
			if isRaDec(line):
				d = spl2[1].split('d') #split out d,m,s
				m = d[1].split('m')
				s = m[1].split('s')
				df,mf,sf = float(d[0]),float(m[0]),float(s[0]) #cast as floats (ints will be recast on assignment later)
				inDeg = df+((mf+(sf/60.))/60.)
				exec(('%s.%s=%f')%(obj,spl2[0],inDeg)) in globals(), locals()
			elif isName(line):
				exec(("%s.%s='%s'")%(obj,spl2[0],spl2[1])) in globals(), locals()
			else:
				exec(("%s.%s=float('%s')")%(obj,spl2[0],spl2[1])) in globals(), locals()
			
		# execute observation assignments
		elif isObservation(line):
			spl = line.split(',') #split out JD,True,Obs
			obsList.append(float(spl[0]))
			RV_True.append(float(spl[1]))
			RV_Obs.append(float(spl[2]))
			
		# update the object we're working with
		elif isNewObject(line):
			if line.startswith('STAR'):
				obj = 'ps' #protostar
			elif line.startswith('REAL PLANET'):
				i = int(line[-1]) - 1
				tryToStore(pp_arr,i,ProtoPlanet())
				obj = 'pp_arr[%d]'%i
			elif line.startswith('FIT PLANET'):
				i = int(line[-1]) - 1
				tryToStore(pp_fit_arr,i,ProtoPlanet())
				obj = 'pp_fit_arr[%d]'%i
	
	#Create Star based on ProtoStar
	coord = SkyCoord(ra=ps.ra*u.degree, dec=ps.dec*u.degree)
	star = Star(ps.massStar,ps.t0,coord,ps.sigma,planets=None,name=ps.name)
	star.nPlanets = int(ps.nPlanets)
	star.nPlanets_fit = int(ps.nPlanets_fit)
	
	#Create Planets based on ProtoPlanets
	assert star.nPlanets == len(pp_arr),'Ensure that you skipped no numbers when numbering the real planets.'
	planets = []
	for pp in pp_arr:
		p = Planet(star,pp.mass,pp.period,pp.ecc,pp.incl,pp.w,pp.tp,K=pp.k,isReal=True)
		planets.append(p)
	star.planets = planets
	star.obsList = obsList
	star.RV_True = RV_True
	star.RV_Obs = RV_Obs
	
	#Create Fitted Planets based on ProtoPlanets
	assert star.nPlanets_fit == len(pp_fit_arr),'Ensure that you skipped no numbers when numbering the fitted planets.'
	planets_fit = []
	for pp in pp_fit_arr:
		p = Planet(star,pp.mass,pp.period,pp.ecc,pp.incl,pp.w,pp.tp,K=pp.k,isReal=False)
		planets_fit.append(p)
	star.planets_fit = planets_fit
	
	return star
					
def printOutput(star):
	# Injected planets' parameters
	inj_str = '\nINJECTED PARAMETERS:'
	for p in star.planets:
		inj_str += '\n'+p.name
		inj_str += ("\nPeriod = {0:.4f} days\nEccentricity = {1:.3f}"
					"\nPeriastron Time = {2:.3f} (days-t0)\nArgument of Periastron = {3:.2f} deg"
					"\nSemiamplitude = {4:.2f} m/s\nMass = {5:.2f} M_e\n").format(*p.params_print)
	# Fitted planets' parameters & Uncertainties
	fit_str = '\nFITTED PARAMETERS:'
	unc_str = ''
	for p in star.planets_fit:
		fit_str += '\n'+p.name
		fit_str += ("\nPeriod = {0:.4f} days\nEccentricity = {1:.3f}"
					"\nPeriastron Time = {2:.3f} (days-t0)\nArgument of Periastron = {3:.2f} deg"
					"\nSemiamplitude = {4:.2f} m/s\nMsini = {5:.2f} M_e\n").format(*p.params_print)
		
		try:
			dum = p.uncs
			for p in star.planets_fit:
				unc_str += '\n'+p.name
				unc_str += ("\n1-sigma uncertainties:\nPeriod = {0:.4f} days\nEccentricity = {1:.3f}"
							"\nPeriastron Time = {2:.3f} (days-t0)\nArgument of Periastron = {3:.2f} rad"
							"\nSemiamplitude = {4:.2f} m/s\nMsini = {5:.2f} M_e\n").format(*p.uncs)
		except AttributeError:
			pass

	print inj_str
	print fit_str
	if unc_str=='':
		print 'Uncertainties not calculated.'
	else:
		print unc_str
	sys.stdout.flush()

def SaveSytem(star,planets):
	return True
