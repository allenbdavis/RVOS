import numpy as np
from tqdm import tqdm
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import AltAz, get_sun
from RVOS_Orbit import calcTrueAnomaly

def observeStar(star,ndays,std_True,t0_in,coord,loc,t_res,objAltThres,weather):

    planets = star.planets
    n_pls = np.size(planets)

    t0 = Time(t0_in,format='jd')
    
    t = t0 #the day considered in each loop
    sep = 1.*u.day
    obsList = []
    
    pbar = tqdm(total=ndays)
    attempted_obs = 0
    
    def pickingLoop(obsList,t,attempted_obs,sep,ndays,flag1=True,flag2=True):
        while (attempted_obs < ndays) and (flag1 or np.size(obsList) < 5) and (flag2 or np.size(obsList) < 10):
            nextObsTime = pickObsTime(coord,loc,t,t_res,objAltThres)
            todaysweather = np.random.random() #how's it looking tonight?
            
            if todaysweather < weather: #it's clear; let's observe!
                obsList = np.append(obsList,nextObsTime)
                pbar.update(1)
            t = nextObsTime + sep
            attempted_obs += 1
        return obsList,t,attempted_obs
    
    # Start with ~1 day separation between observations
    obsList,t,attempted_obs = pickingLoop(obsList,t,attempted_obs,sep,ndays,flag1=False)
    
    # Now decide whether to use a longer spacing
    if scatter(star,obsList) > std_True: # Keep a separation of 1 day
        obsList,t,attempted_obs = pickingLoop(obsList,t,attempted_obs,sep,ndays) 

    else: # Ramp it up to 5 day separation
        sep = 5.*u.day
        obsList,t,attempted_obs = pickingLoop(obsList,t,attempted_obs,sep,ndays,flag2=False)  
        
        # Reevaluate whether the separation must again be increased
        if scatter(star,obsList) > std_True: # Keep a 5 day separation
            obsList,t,attempted_obs = pickingLoop(obsList,t,attempted_obs,sep,ndays)
            
        else: # Ramp it up to 10 day separation and keep it there
            sep = 10.*u.day
            obsList,t,attempted_obs = pickingLoop(obsList,t,attempted_obs,sep,ndays)

    pbar.close()
    # Save the obsList
    star.obsList_units = obsList #in terms of astropy units
    star.obsList = removeUnits(obsList) #simple JD values

    # Record true RV curve
    RVs = [None]*n_pls
    for i in range(0,n_pls):
        p = planets[i]
        f = calcTrueAnomaly(p.period,p.tp,p.ecc,removeUnits(obsList))
        RVs[i] = p.K*(np.cos(p.w+f) + p.ecc*np.cos(p.w))
        p.RV_True = RVs[i] #also save each individual planet's RV
        
    RV = np.sum(RVs,axis=0) #sum each planet's RV contribution
    star.RV_True = RV
    
    #add jitter scaled to a normal distribution
    errList = np.random.normal(0, std_True, np.size(star.obsList))

    #record observed RV curve
    star.RV_Obs = star.RV_True + errList

    #record observated scatter in RV curve
    star.std_Obs = np.std(star.RV_Obs)
    


# Calculates all the possible times a star will be above an altitude and returns a random time
# to observe it while the sun is down
def pickObsTime(coord,loc,t,t_res,objAltThres,sunAltThres=-18.):
    #coord is SkyCoord object; loc is EarthLocation object; t0 is Time object
    #t_res is time resolution in minutes. Suggest 5 or 10 minutes
    
    delta_t = np.linspace(0,24,24*60./t_res)*u.hour #divides a day into 5-minute chunks
    timeRange = t+delta_t
        
    altazframe = AltAz(obstime=timeRange, location=loc)
    sunaltazs = get_sun(timeRange).transform_to(altazframe)
    objaltazs = coord.transform_to(altazframe)
    
    sunDownTimes = delta_t[np.where(sunaltazs.alt < sunAltThres*u.deg)]
    objUpTimes = delta_t[np.where(objaltazs.alt > objAltThres*u.deg)]

    sunBools = sunaltazs.alt < sunAltThres*u.deg #True/False arrays
    objBools = objaltazs.alt > objAltThres*u.deg
    
    goodoffsets = delta_t[sunBools*objBools]
    size = np.size(goodoffsets)
    
    if size == 0:
        return pickObsTime(coord,loc,findNextOpportunity(t,coord,loc,t_res,objAltThres),
                           t_res,objAltThres)
    else:
        goodtimes = t+goodoffsets #picks out the times that satisfy both conditions
        return goodtimes[np.random.randint(0,size)]

    
#returns the date when the object will next be good to observe
def findNextOpportunity(t,coord,loc,t_res,objAltThres):
    start = t
    end = start + 1.*u.year
    found = False
    while start <= end and not found:
        midpoint = Time((start.value+end.value)/2,format='jd')
        if isGood(midpoint,coord,loc,t_res,objAltThres): 
            end = midpoint
            if not isGood(midpoint-1*u.day,coord,loc,t_res,objAltThres):
                found = True
                return midpoint
        else:
            start = midpoint + 1*u.day
    return midpoint


# tests whether the object is good to observe on a certain JD
def isGood(t,coord,loc,t_res,objAltThres,sunAltThres=-18.):
    delta_t = np.linspace(0,24,24*60./t_res)*u.hour #divides a day into 5-minute chunks
    timeRange = t+delta_t

    altazframe = AltAz(obstime=timeRange, location=loc)
    sunaltazs = get_sun(timeRange).transform_to(altazframe)
    objaltazs = coord.transform_to(altazframe)

    sunDownTimes = delta_t[np.where(sunaltazs.alt < sunAltThres*u.deg)]
    objUpTimes = delta_t[np.where(objaltazs.alt > objAltThres*u.deg)]

    sunBools = sunaltazs.alt < sunAltThres*u.deg #True/False arrays
    objBools = objaltazs.alt > objAltThres*u.deg

    goodoffsets = delta_t[sunBools*objBools]
    size = np.size(goodoffsets)

    return not size == 0

def scatter(star,obsList):
    planets = star.planets
    n_pls = np.size(planets)
    
    RVs = [None]*n_pls
    #find the RV contribution of each planet independently
    for i in range(0,n_pls):
        p = planets[i]
        f = calcTrueAnomaly(p.period,p.tp,p.ecc,removeUnits(obsList))
        RVs[i] = p.K*(np.cos(p.w+f) + p.ecc*np.cos(p.w))
    
    RV = np.sum(RVs,axis=0) #sum each planet's RV contribution
    return np.std(RV)

def removeUnits(arr_in):
    if not isinstance(arr_in[0],float):
        arr = np.array(arr_in[0].value)
        for elem in arr_in[1:]:
            arr = np.append(arr,elem.value)
    return arr