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
#




def main():

    start_time = datetime.now()
        
    ###---  differential problem  ---###
    eq1 = lambda t,u : u[1]
    eq2 = lambda t,u : u[2]
    eq3 = lambda t,u : 4*u[0] - 3 *u[1] + 2*u[2]

    func1 = np.array([eq1,eq2,eq3])
    
    y0      = np.array([0.1,0.5,0.1])
    system1 = rhs.Rhs(func1, 0,2,y0,100 )

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
#----------------------------------------------------------------------------------------------   
    end_time = datetime.now()
    print('Duration {}'.format(end_time - start_time))
    
    fig,axs = PalatinoItalic()(1,1,figsize=(9.5,4.5))
    axs.plot(fet,feu,label='Exp Euler')
    axs.plot(iet,ieu,label='Imp Euler')
    #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
    #axs.plot(bet,beu[:,3],label='NR-BWDEuler') #, linestyle=':')
    #axs.plot(lft,lfu,label='Leap Frog')
    #axs.plot(cnt,cnu,label='Crank Nicholson')
    #axs.plot(amt,amu,label='Adams Moulton 5th') #,linestyle=':')
    #axs.plot(rkt,rku,label='Runge Kutta 4th',linestyle=':')

    legend = leg = axs.legend()
    legend.get_frame().set_linewidth(1.2)
    plt.setp(legend.get_texts(), color='#555555')
    axs.set(xlabel = 'time' ,ylabel='y(t)')
    axs.set_title('dy$_1$/dt = y$_2$  ;  dy$_2$/dt = y$_3$ ; dy$_3$/dt = 4y$_1$+3y$_2$+2y$_3$  ',y=1.05)
    plt.savefig('ODE/Results/system04.pdf')   
    plt.show()


if __name__ == '__main__':
    main()









