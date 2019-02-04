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
    
    eq1 = lambda t,u : (2-0.5*u[1])*u[0]   
    eq2 = lambda t,u : (-1+0.5*u[0])*u[1]   
        
    func1 = np.array([eq1,eq2])
    
    y0      = np.array([6.,2.])
    system1 = rhs.Rhs(func1, 0,15,y0,1200 )

    fwdeuler = euler.Explicit(system1)
    fet,feu  = fwdeuler.solve()   

    bwdeuler = euler.Backward(system1)
    bet,beu  = bwdeuler.solve()
    
    sieuler = euler.SemiImplicit(system1)
    siet,sieu  = sieuler.solve()
    
    Impeuler = euler.Implicit(system1)
    iet,ieu  = Impeuler.solve()
    

    #leapfrog    = multistep.LeapFrog(system1, 'lf1.dat')
    #lft,lfu   = leapfrog.solve()
    
#    cranknicholson = rungekutta.CrankNicholson(problem1 , 'cn1.dat')
#    cnt,cnu    = cranknicholson.solve()
#    
#    am5 = adamsmethods.AdamsMoulton(problem1, 'AM5_1.dat')
#    amt,amu = am5._5th.solve()
#        
    rk3 = rungekutta.RK3(system1,'s1RK4.dat')
    rkt,rku = rk3.solve()
#-----------------------------------------------------------------------------------------------------   
    end_time = datetime.now()
    print('Duration {}'.format(end_time - start_time))
    #**{'color': 'paraview'}
    fig,axs = Standard()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    
    #fig,axs = qualityplot.StandardImproved(figsize=(9.5,4.5),**{'scheme':'ernest'})(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Elsevier(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexMathpazo(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexMathpazoBeamer(figsize=(6.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TexPaper(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Beamer(figsize=(6.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Helvet(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Times(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TimesItalicDefault(figsize=(9.5,4.5))(1,1,**{ 'legend.loc':'upper left'})
    #fig,axs = qualityplot.TimesItalicImproved(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.TimesItalicModified(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.SansModified(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.Fonts(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = qualityplot.PalatinoItalic(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    
    #axs.plot(x,y,label='Analitycal', marker='o',linestyle='',markersize=4.5,color='C0',alpha=0.7)
    axs.plot(fet,feu,label='Exp Euler')
    #axs.plot(iet,ieu,label='Imp Euler')
    #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
    axs.plot(bet,beu,label='NR-BWDEuler') #, linestyle=':')
    #axs.plot(lft,lfu,label='Leap Frog')
    #axs.plot(cnt,cnu,label='Crank Nicholson')
    #axs.plot(amt,amu,label='Adams Moulton 5th') #,linestyle=':')
    axs.plot(rkt,rku,label='Runge Kutta 3th')

    legend = leg = axs.legend()
    legend.get_frame().set_linewidth(1.2)
    plt.setp(legend.get_texts(), color='#555555')
    axs.set(xlabel = 'time' ,ylabel='y(t)')
    axs.set_title('dy$_1$/dt= -4 (y$_1$ - y$_1$y$_2$)  ;  dy$_2$/dt = -1(y$_2$ - y$_1$y$_2$)',y=1.05)
    plt.savefig('ODE/Results/system06b.pdf')   
    plt.show()


if __name__ == '__main__':
    main()









