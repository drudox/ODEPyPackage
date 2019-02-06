#!/usr/bin/env python3
# -*- coding : utf-8 -*-



import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import sys

sys.path.append('../../')

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
   

    eq1 = lambda t,u : u[1]*u[2] 
    eq2 = lambda t,u : -2*u[0]*u[2]    
    eq3 = lambda t,u : u[0]*u[1]

    func1 = np.array([eq1,eq2,eq3])
    
    y0      = np.array([0.1039,0.9946,-0.0069])
    system1 = rhs.Rhs(func1, 0,20,y0,100 )

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
    #fig,axs = Standard(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = StandardImproved(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = Elsevier(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TexMathpazo(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TexMathpazoBeamer(figsize=(6.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TexPaper(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = Beamer(figsize=(6.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = Helvet(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = Times(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TimesItalicDefault(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TimesItalicImproved(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TimesItalicModified(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    fig,axs = SansItalic()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})
    #fig,axs = Fonts(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = PalatinoItalic(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    
    #axs.plot(x,y,label='Analitycal', marker='o',linestyle='',markersize=4.5,color='C0',alpha=0.7)
    axs.plot(fet,feu,label='Exp Euler')
    #axs.plot(iet,ieu,label='Imp Euler')
    #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
    axs.plot(bet,beu,label='NR-BWDEuler') #, linestyle=':')
    #axs.plot(lft,lfu,label='Leap Frog')
    #axs.plot(cnt,cnu,label='Crank Nicholson')
    #axs.plot(amt,amu,label='Adams Moulton 5th') #,linestyle=':')
    #axs.plot(rkt,rku,label='Runge Kutta 3th')

    legend = leg = axs.legend()
    legend.get_frame().set_linewidth(1.2)
    plt.setp(legend.get_texts(), color='#555555')
    axs.set(xlabel = 'time' ,ylabel='y(t)')
    axs.set_title('dy$_1$/dt= y$_2$ y$_3$ ; dy$_2$/dt= -2y$_1$y$_3$ ; dy$_3$/dt = y$_1$y$_2$',y=1.05)
    plt.ylim(-1.5,1.7)
    plt.savefig('../Results/system09.pdf')   
    plt.show()


if __name__ == '__main__':
    main()









