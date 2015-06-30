def SaveSystem(star,planets,fname=''):
	from time import strftime
	import os
	import errno
	s=''
	s+='#This file contains the info necessary to recreate the 	Star, Planets, and observation data.\n'
	s+='\nSTAR\n'
	s+='name=%s\n'%star.name
	s+='massStar=%s\n'%star.massStar
	s+='t0=%s\n'%star.t0
	s+='RA=%s\n'%star.coord.ra
	s+='DEC=%s\n'%star.coord.dec
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
		s+='K=%s\n'%pl.K
		s+='isReal=%s\n'%pl.isReal
	
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
		s+='K=%s\n'%pl.K
		s+='isReal=%s\n'%pl.isReal
		
		try:
			dum = pl.uncs
			s+='UNCS EXIST\n'
			s+='uncs[5]=%s\n'%pl.uncs[5] #mass
			s+='uncs[0]=%s\n'%pl.uncs[0] #period
			s+='uncs[1]=%s\n'%pl.uncs[1] #ecc
			s+='uncs[3]=%s\n'%pl.uncs[3] #w
			s+='uncs[2]=%s\n'%pl.uncs[2] #tp
			s+='uncs[4]=%s\n'%pl.uncs[4] #K
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

def LoadSystem(fname):
	f = open(fname,'r')
	s = f.readlines()
	f.close()
	
	star = Star
	
	def isAssignment(line):
    #lines start with a lower case letter IFF they are assignment lines
		return line[0].islower()

	def isNewObject(line):
		#lines start with an upper case letter IFF they are new object lienes
		return line[0].isupper()
		
	def isObservation(line):
		#lines start with a numeral IFF they are observation data
		return line[0].isdigit()

	
	obj = None
	star.obsList = []
	star.RV_True = []
	star.RV_Obs = []
	for line in s:
		if not obsMode:
			# execute a parameter assignment
			if isAssignment(line):
				spl1 = line.split('#') #split out comments
				spl2 = spl1[0].split('=') #split apart assignment statement
				assert not obj==None,'OBJECT must be specified before its parameters.'
				exec(('%s.%s=%s')%(obj,spl2[0],spl2[1])) #execute the assignment
			
			# execute observation assignments
			elif isObservation(line):
				spl = line.split(',') #split out JD,True,Obs
				exec('star.obsList.append(%s)'%spl[0]) #execute the assignments
				exec('star.RV_True.append(%s)'%spl[1])
				exec('star.RV_Obs.append(%s)'%spl[2])
				
			# update the object we're working with
			elif isNewObject(line):
				if line=='STAR':
					obj = 'star'
				elif line.startswith('REAL PLANET'):
					i = line[-1] - 1
					obj = 'star.planets[%d]'%i
				elif line.startswith('FITTED PLANET'):
					i = line[-1] - 1
					obj = 'star.planets_fit[%d]'%i