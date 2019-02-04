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
    
          
    # differential problem 
    
    a = 0.1
    b = 0.001
    c = 0.0012
    d = 0.08

    eq1 = lambda t,u : u[0] * (a-b*u[1]) 
    eq2 = lambda t,u : u[1] * (d*u[0] - c) 
        
    func1 = np.array([eq1,eq2])
    
    y0      = np.array([10.,0.4])
    system1 = rhs.Rhs(func1, 0,10,y0,100 )

    fwdeuler = euler.Explicit(system1)
    fet,feu  = fwdeuler.solve()   

    bwdeuler = euler.Backward(system1)
    bet,beu  = bwdeuler.solve()
    
    sieuler = euler.SemiImplicit(system1)
    siet,sieu  = sieuler.solve()
    
   # Impeuler = euler.Implicit(system1, 's1bwe1.dat')
   # iet,ieu  = Impeuler.solve()
    

    #leapfrog    = multistep.LeapFrog(system1, 'lf1.dat')
    #lft,lfu   = leapfrog.solve()
    
#    cranknicholson = rungekutta.CrankNicholson(problem1 , 'cn1.dat')
#    cnt,cnu    = cranknicholson.solve()
#    
#    am5 = adamsmethods.AdamsMoulton(problem1, 'AM5_1.dat')
#    amt,amu = am5._5th.solve()
#        
    rk3 = rungekutta.RK3(system1)
    rkt,rku = rk3.solve()

#-----------------------------------------------------------------------------------------------------   
    end_time = datetime.now()
    print('Duration {}'.format(end_time - start_time))
            #**{'color': 'paraview'}
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
    fig,axs = PalatinoItalic()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})

   






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
    axs.set_title('dy$_1$/dt= y$_1$ (a-b y$_2$) ;  dy$_2$/dt =y$_2$ (d y$_1$ -c)' ,y=1.05)
    plt.ylim(-50,1700)
    plt.savefig('ODE/Results/system08.pdf')   
    plt.show()


if __name__ == '__main__':
    main()









