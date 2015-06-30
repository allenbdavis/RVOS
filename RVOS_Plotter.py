import matplotlib.pyplot as plt
import numpy as np
import os
import errno
from RVOS_Orbit import calcTrueAnomaly
from time import strftime

def makePlot(sim,numpoints=500,fill=True,save='',fmt='png',hide_True=False,date=False):
    
    star = sim.star
        
    if fill:
        # get more data points for nicer plotting
        start,end = star.obsList[0],star.obsList[-1]
        
        if numpoints == -1:
            # default is 10 points per period of fastest planet
            # find fastest period planet, among the fitted and real planets
            fastestperiod = 10000000
            for p in star.planets_fit:
                if p.period < fastestperiod: fastestperiod = p.period
            for p in star.planets:
                if p.period < fastestperiod: fasterperiod = p.period
            numpoints = int(np.round(10*(end-start)/fastestperiod))
            if numpoints > 3000: numpoints = 3000

        times = np.linspace(start,end,numpoints)

        if not hide_True:
            RV_True = [None]*star.nPlanets
            for i in range(0,star.nPlanets):
                p = star.planets[i]
                f_True = calcTrueAnomaly(p.period,p.tp,p.ecc,times)
                RV_True[i] = p.K*(np.cos(p.w+f_True) + p.ecc*np.cos(p.w))
            RV_True = np.sum(RV_True,axis=0)

        RV_fit = [None]*star.nPlanets_fit
        for i in range(0,star.nPlanets_fit):
            p = star.planets_fit[i]
            f_fit = calcTrueAnomaly(p.period,p.tp,p.ecc,times)
            RV_fit[i] = p.K*(np.cos(p.w+f_fit) + p.ecc*np.cos(p.w))
        RV_fit = np.sum(RV_fit,axis=0)
    
    # Figure parameters
    plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    
    #Plot the RV curve
    plt.subplot(2, 1, 1)

    plt.plot(star.obsList,star.RV_Obs,'kx')
    
    if fill:
        plt.plot(times, RV_fit, 'b--',label='Fit')
        if not hide_True: plt.plot(times, RV_True,'r--',label='True')
    else:
        plt.plot(star.obsList, star.RV_fit,'b-',label='Fit')
        if not hide_True: plt.plot(star.obsList, star.RV_True,'r--',label='True')

#     plt.errorbar(star.obsList, star.RV_Obs, yerr=star.std_True, capsize=0,linestyle = 'none',ecolor='k')
    plt.xlabel('Time (days)')
    plt.ylabel('RV (m/s)')
    plt.xlim(star.obsList[0],star.obsList[-1])    
    
    
    if sim.find_uncs:
        #ALSO A PAIN WITH MULTIPLE PLANETS. DISABLING FOR NOW
        title_str = (r'$M_{star}$=%.2f $M_{\odot}$, $N_{obs}$=%d, $\sigma_{RV}$=%.2f m/s'
                     '\nTrue: per=%.4f days, mass=%.2f $M_{\oplus}$, ecc=%.2f, w=%.2f deg'
                     '\nFit: per=%.4f$\pm$%.4f days, msini=%.2f$\pm$%.2f $M_{\oplus}$,'
                     'ecc=%.2f$\pm$%.2f, w=%.2f$\pm$%.2f deg'
                     '')%(p.massStar,np.size(p.obsList),star.sigma,
                          p.period,p.mass,p.ecc,p.w*180./np.pi,
                          p.params_out_print[0],p.uncs[0],p.params_out_print[5],p.uncs[5],
                          p.params_out_print[1],p.uncs[1],p.params_out_print[3],p.uncs[3]*180./np.pi)
    else:
        title_str = (r'$M_{star}$=%.2f $M_{\odot}$, $N_{obs}$=%d, $\sigma_{RV}$=%.2f m/s, $n_{pl}$=%d'
                     '\n')%(star.massStar,np.size(star.obsList),star.sigma,star.nPlanets)
        if not star.name == '':
            title_str = star.name+'\n'+title_str
        title_str += 'Real Planets:\n'
        for i in range(0,star.nPlanets):
            p = star.planets[i]
            title_str += (r'%s: P=%.4f days, m=%.2f $M_{\oplus}$, e=%.2f ,w=%.2f$^{\circ}$'
                         '\n')%(chr(i+ord('b')),p.period,p.mass,p.ecc,p.w*180./np.pi)
        title_str += 'Fitted Planets:\n'
        for i in range(0,star.nPlanets_fit):
            p = star.planets_fit[i]
            title_str += (r'%s: P=%.4f days, msini=%.2f $M_{\oplus}$, e=%.2f, w=%.2f$^{\circ}$'
                         '\n')%(chr(i+ord('b')),p.period,p.mass,p.ecc,p.w*180./np.pi)   
    
    plt.title(title_str,fontsize=10)
    plt.legend(fontsize=10)
    
    #Plot the Periodogram
    plt.subplot(2, 1, 2)

    plt.plot(sim.periodList, star.P_G, 'k-')
    
    # For each fitted planet
    for p in star.planets_fit:
        plt.axvline(p.period,ymin=.95,color='b',linewidth=.5)
    # For each real planet
    for p in star.planets:
        plt.axvline(p.period,ymin=.95,color='r',linewidth=.5)
    
    # plot the significance lines
    if sim.find_FAP:
        xlim = (np.min(sim.periodList), np.max(sim.periodList))
        for pair in star.FAP_powers:
            lvl,power = pair
            plt.plot(xlim,[power,power],'--', lw=1., label='%.2f%%'%lvl)
    
#     plt.xlim(sim.freqList[0],sim.freqList[-1])
    
    plt.xlabel('Period (days)')
    plt.ylabel('Power')
    plt.xscale('log')
#     plt.title('Ntrials = %d' %sim.FAP_niter, fontsize=12)
    plt.tight_layout()
    plt.legend(fontsize=11,loc=2)
    
    
    
    if save != '':
        if hide_True:
            save += '_hide'
        
        # Make a new directory every day
        path = '.\\figures_'+strftime("%Y-%m-%d")
        
        
        def make_sure_path_exists(path):
            try:
                os.makedirs(path)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            
        make_sure_path_exists(path)
        fmts = fmt.split(',')
        for f in fmts:
            plt.savefig(path+'\\'+save+'_'+sim.theTime+'.'+f)
    else:
        plt.show()
		
def makePhasePlot(star):
    
    #Plot a phase folded version
    
    assert star.nPlanets == star.nPlanets_fit, 'makePhasePlot requires you to fit as many planets as there actually are.'
    
    def phaseFold(times, period):
        return (times % period)/period

    for i in range(0,star.nPlanets):
        p = star.planets[i]
        p_fit = star.planets_fit[i]
        plt.subplot(star.nPlanets, 1, i)
        plt.plot(phaseFold(star.obsList, p.period),star.RV_True,'r.',label='True')
        plt.plot(phaseFold(star.obsList, p_fit.period),star.RV_Obs,'b.',label='Observed')
        plt.xlabel('Phase (given Period = %.4f days)' %p.period )
        plt.ylabel('RV (m/s)')
        if not p.name == '':
            plt.title(p.name)
        plt.tight_layout()
        plt.legend()
    
    
#     plt.subplot(2,1,2)
#     for p in star.planets:
#     plt.plot(phaseFold(star.obsList, p.params_out[0]),p.RV_True,'r.',label='True')
#     plt.plot(phaseFold(star.obsList, p.params_out[0]),p.RV_Obs,'b.',label='Observed')
#     plt.xlabel('Phase (given Period = %.4f days)' %p.params_out[0] )
#     plt.ylabel('RV (m/s)')
#     plt.title('Fitted period')
#     plt.tight_layout()
#     plt.legend()

    plt.show()