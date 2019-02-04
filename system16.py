#!/usr/bin/env python3
# -*- coding : utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

from ODE.Rhs import rhs
from ODE.Solvers.Euler import euler
from ODE.Solvers.RungeKutta import rungekutta
from ODE.Solvers.MultiStep import multistep
from ODE.Solvers.AdamsMethods import adamsmethods
from ODE.qualityPlot.qualityplot import *
from ODE.Solvers.CentDiff import centdiff
#


def main():

    start_time = datetime.now()
    #############################################       
    # Chemostats and Microorganism Culturing 
    #############################################
    eq11 = lambda t,u : -0.075*u[0] - 0.008177008175*u[1]     
    eq21 = lambda t,u :  0.4401515152 *u[1] 
    
    eq12 = lambda t,u: -1.690372243*u[0] - 0.001190476190 * u[1]
    eq22 = lambda t,u: 101.7684513*u[0] 

    func1 = np.array([eq11,eq21])
    func2 = np.array([eq12,eq22])


    y01      = np.array([20.5,0.4])
    y02      = np.array([0.3,11.4])
    
    system1 = rhs.Rhs(func1,0,10,y01,500 )
    system2 = rhs.Rhs(func2,0,10,y02,500 )

    fwdeuler = euler.Explicit(system1)
    fe1t,fe1u  = fwdeuler.solve()   

    bwdeuler = euler.Backward(system1)
    be1t,be1u  = bwdeuler.solve()
  
    fwdeuler = euler.Explicit(system2)
    fe2t,fe2u  = fwdeuler.solve()   

    bwdeuler = euler.Backward(system2)
    be2t,be2u  = bwdeuler.solve()
    
#    sieuler = euler.SemiImplicit(system2)
#    siet,sieu  = sieuler.solve()
    
    #Impeuler = euler.Implicit(system1)
    #iet,ieu  = Impeuler.solve()
    
#    midpoint = centdiff.MidPoint(system1)
#    cdt,cdu  = midpoint.solve()
    
    #leapfrog    = multistep.LeapFrog(system1)
    #lft,lfu   = leapfrog.solve()
    
#    cranknicholson = rungekutta.CrankNicholson(system1)
#    cnt,cnu    = cranknicholson.solve()
#    
#    am5 = adamsmethods.AdamsMoulton(system1)
#    amt,amu = am5._5th.solve()
#        
#    rk4 = rungekutta.RK4(system1)
#    rkt,rku = rk4.solve()

#-----------------------------------------------------------------------------------------------------   
    end_time = datetime.now()
    print('Duration {}'.format(end_time - start_time))
            #**{'color': 'paraview'}
    #fig,axs = qualityplot.Standard(**{'scheme':'vega'})(1,1,figsize=(8.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.StandardImproved(**{'scheme':'vega'} )(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Elsevier()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Elsevier(figsize=(9.5,4.5),**{'scheme':'mycolor'})(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexMathpazo(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexMathpazo(figsize=(9.5,4.5),**{'scheme':'mycolor'})(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexMathpazoBeamer(figsize=(6.5,4.5),**{'scheme':'mycolor'})(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexPaper()(1,1,figsize=(8.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Helvet()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.HelvetBeamer()(1,1,figsize=(8.5,5.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Beamer()(1,1,figsize=(8.5,5.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Times()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TimesItalicDefault(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TimesItalicImproved()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TimesItalicModified()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.SansModified(figsize=(9.5,4.5),**{'scheme':'vega'} )(1,1)
    #fig,axs = qualityplot.Fonts(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Palatino(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    fig,axs = PalatinoItalic()(1,1, figsize=(9.5,4.5) )# ,**{'scheme':'nb'})
    
    #axs.plot(x,y,label='Analitycal', marker='o',linestyle='',markersize=4.5,color='C0',alpha=0.7)
    axs.plot(fe1t,fe1u,label='Exp Euler 1')#,color='C0')
    #axs.plot(iet,ieu,label='Imp Euler')
    #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
    #axs.plot(be1t,be1u,label='NR-BWDEuler 1', color='C1') #, linestyle=':')
    #axs.plot(lft,lfu,label='Leap Frog')
    #axs.plot(cnt,cnu,label='Crank Nicholson')
    #axs.plot(amt,amu,label='Adams Moulton 5th') #,linestyle=':')
    #axs.plot(be1t,be1u,label='NR-BWDEuler 2',color='C1')
    axs.plot(fe2t,fe2u,label='Exp Euler 2')#, color='C1')
    
    chartBox = axs.get_position()    
    legend = leg = axs.legend()
    legend.get_frame().set_linewidth(0.7)
    axs.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
    plt.setp(legend.get_texts(), color='#2E3436')
    plt.legend(bbox_to_anchor=(1.03, 0.75), loc=2, borderaxespad=0.)
    


    #legend = leg = axs.legend()
    #legend.get_frame().set_linewidth(0.7)
    #plt.setp(legend.get_texts(), color='#2E3436')
    axs.set(xlabel = 'time' ,ylabel='y(t)')
    axs.set_title('Chemostats and Micro-organism Culturing ',y=1.05)
    plt.savefig('ODE/Results/system16.pdf')   
    plt.show()


if __name__ == '__main__':
    main()









