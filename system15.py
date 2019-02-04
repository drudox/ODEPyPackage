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
    # differential problem Pond Pollution 
    #############################################
    eq1 = lambda t,u : 0.001*u[2] - 0.001*u[0] + 0.125       
    eq2 = lambda t,u : 0.001*u[0] - 0.001 *u[1]    
    eq3 = lambda t,u : 0.001*u[1] - 0.001 *u[2]     
        
    func1 = np.array([eq1,eq2,eq3])
    
    y0      = np.array([-200.,-140.,10.])
    system1 = rhs.Rhs(func1, 0,2880,y0,500 )

    fwdeuler = euler.Explicit(system1)
    fet,feu  = fwdeuler.solve()   

    bwdeuler = euler.Backward(system1)
    bet,beu  = bwdeuler.solve()
    
    sieuler = euler.SemiImplicit(system1)
    siet,sieu  = sieuler.solve()
    
    Impeuler = euler.Implicit(system1)
    iet,ieu  = Impeuler.solve()
    
    midpoint = centdiff.MidPoint(system1)
    cdt,cdu  = midpoint.solve()
    
    #leapfrog    = multistep.LeapFrog(system1)
    #lft,lfu   = leapfrog.solve()
    
#    cranknicholson = rungekutta.CrankNicholson(system1)
#    cnt,cnu    = cranknicholson.solve()
#    
#    am5 = adamsmethods.AdamsMoulton(system1)
#    amt,amu = am5._5th.solve()
#        
    rk3 = rungekutta.RK3(system1)
    rkt,rku = rk3.solve()

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
    fig,axs = PalatinoItalic()(1,1, figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    
    
    #axs.plot(x,y,label='Analitycal', marker='o',linestyle='',markersize=4.5,color='C0',alpha=0.7)
    axs.plot(fet,feu,label='Exp Euler')
    #axs.plot(iet,ieu,label='Imp Euler')
    #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
    axs.plot(bet,beu,label='NR-BWDEuler') #, linestyle=':')
    #axs.plot(lft,lfu,label='Leap Frog')
    #axs.plot(cnt,cnu,label='Crank Nicholson')
    #axs.plot(amt,amu,label='Adams Moulton 5th') #,linestyle=':')
    #axs.plot(rkt,rku,label='Runge Kutta 3th')
    axs.plot(cdt,cdu,label='CentDiff 2nd')
    
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
    axs.set_title('Pond Pollution',y=1.05)
    plt.savefig('ODE/Results/system15.pdf')   
    plt.show()


if __name__ == '__main__':
    main()









