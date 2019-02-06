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
#



def main():

    start_time = datetime.now()
    a=4;
    c=1;
        
    # differential problem 
    eq1 = lambda t,u : a*(u[0]-u[0]*u[1]);
    eq2 = lambda t,u : -c*(u[1]-u[0]*u[1]);
        
    func1 = np.array([eq1,eq2])
    
    y0      = np.array([2.,1.])
    system1 = rhs.Rhs(func1, 0,10,y0,500 )

    fwdeuler = euler.Explicit(system1)
    fet,feu  = fwdeuler.solve()   

    bwdeuler = euler.Backward(system1)
    bet,beu  = bwdeuler.solve()
    
    sieuler = euler.SemiImplicit(system1)
    siet,sieu  = sieuler.solve()
    
    Impeuler = euler.Implicit(system1)
    iet,ieu  = Impeuler.solve()
    

    leapfrog    = multistep.LeapFrog(system1)
    lft,lfu   = leapfrog.solve()
        
    rk3 = rungekutta.RK3(system1)
    rkt,rku = rk3.solve()

    fig,axs = SansItalic()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})


#----------------------------------------------------------------------------------------------   
    end_time = datetime.now()
    print('Duration {}'.format(end_time - start_time))
            #**{'color': 'paraview'}
    #fig,axs = Standard(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = StandardImproved(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = Elsevier(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TexMathpazo(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TexMathpazoBeamer(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = Beamer(figsize=(6.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TexPaper(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = Helvet(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = HelvetBeamer(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = Times(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TimesItalicDefault(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = TimesItalicImproved(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = PalatinoItalic(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})
    #fig,axs = Fonts(figsize=(9.5,4.5))(1,1)# ,**{'scheme':'nb'})

    axs.plot(fet,feu,label='Exp Euler')
    #axs.plot(iet,ieu,label='Imp Euler')
    #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
    axs.plot(bet,beu,label='NR-BWDEuler') #, linestyle=':')
    axs.plot(lft,lfu,label='Leap Frog')
    #axs.plot(cnt,cnu,label='Crank Nicholson')
    #axs.plot(amt,amu,label='Adams Moulton 5th') #,linestyle=':')
    #axs.plot(rkt,rku,label='Runge Kutta 3th',linestyle=':')

    legend = leg = axs.legend()
    legend.get_frame().set_linewidth(1.2)
    plt.legend(loc='upper left')
    #axs.plot(x,y,label='Analitycal', marker='o',linestyle='',markersize=4.5,color='C0',alpha=0.7)
    plt.setp(legend.get_texts(), color='#555555')
    axs.set(xlabel = 'time' ,ylabel='y(t)')
    axs.set_title('dy$_1$/dt= -4 (y$_1$ - y$_1$y$_2$)  ;  dy$_2$/dt = -1(y$_2$ - y$_1$y$_2$)',y=1.05)
    plt.savefig('../Results/system02.pdf')   
    plt.show()


if __name__ == '__main__':
    main()









