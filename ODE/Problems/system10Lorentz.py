#!/usr/bin/env python3

'''   
  Driver programs for testing solvers for: 
  Sistem of first order differential equations 
  
  @author Marco Ghiani Glasgow 2018

  developing for python v.3.x 
'''



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




#------------------------------------------------------------
## Lorentz attractor system

func1a = lambda t,u : 10   * (u[1] - u[0])
func2a = lambda t,u : 28   * u[0] - u[1] - u[0] * u[2]  
func3a = lambda t,u : -8/3 * u[2] + u[0]*u[1]

#-------------------------------------------------------------
#-------------------------------------------------------------
func1 = np.array([func1a,func2a,func3a])
#

def main():
   
   y01 = np.array([1.,0.,0.])
   
   problem1 = rhs.Rhs(func1, 0.0, 50.0 , y01, 10000)
   
   #lp_p1    = multistep.LeapFrog(problem1)
   #lft,lfu  = lp_p1.solve()
#    
#
#
#  ab2_p1 = adamsmethods.AdamsBashforth(problem1)
   #ab2_p  = ab2_p1._2nd(problem1)
#  ab2t,ab2u = ab2_p1._2nd.solve()
   #ab2_p1._2nd.plot()
 
#   am2_p1 = adamsmethods.AdamsMoulton(problem1)
   #ab2_p  = ab2_p1._2nd(problem1)
#   am2t,am2u = am2_p1._2nd.solve()
   #ab2_p1._2nd.plot()


   #rkt =
#   rkt = np.genfromtxt('rk2_1.dat',usecols=(0,)) 
#   rku = np.genfromtxt('rk2_1.dat',usecols=(1,)) 
 

   #print(type(rk_t) , type(rk_u), type(cnu) ) 
   #print(rk_t.shape , rk_u.shape, cnt.shape ) 
#   
#   bwdeuler_p4 = euler.Implicit(problem4, 'bwe.dat')
#   bet4,beu4   = bwdeuler_p4.solve()
#
#   cn_p2      = rungekutta.CrankNicholson(problem1, 'CN_lorentz.dat')
#   cnt,cnu    = cn_p2.solve()
#    
#   rk_p2      =  rungekutta.RK2.Heun(problem1, 'Heun_lorentz.dat')
#   rket,rkeu  =  rk_p2.solve()
#
#   fwdeuler_p1 = euler.Explicit(problem1, 'FWE_lorentz.dat')
#   fet,feu = fwdeuler_p1.solve()

   #bwdeuler_p1 = euler.Implicit(problem1, 'BWE_lorentz.dat')
   #bet,beu = bwdeuler_p1.solve()
   
   #bwdeuler_p1 = euler.Backward(problem1, 'BWE_lorentz.dat')
   #bet,beu = bwdeuler_p1.solve()

   rk4 = rungekutta.RK4(problem1)
   rkt,rku = rk4.solve()
 

#   #leapfrog_p1 = euler.Explicit(problem1, True, True, 'lorentz_lf.dat')
#   #lf1,lf1     = leapfrog_p1.solve()
#  
#   leapfrog_p2 = euler.Explicit(problem2, 'oscillator_lf.dat')
#   lf2,lf2     = leapfrog_p2.solve()





   #fwdeuler_p1.plot()
   
   #ab2_p1 = euler.Implicit(problem3,True,True,filename='lorentz_bwe.dat')
   #bet,beu = bwdeuler_p1.solve()
   #bwdeuler_p1.plot()
   
   
   

   fig,axs = PalatinoItalic()(1,1,figsize=(9.5,4.5))# ,**{'scheme':'nb'})

   #axs.plot(fet,feu,label='Exp Euler')
   #axs.plot(iet,ieu,label='Imp Euler')
   #axs.plot(bet,beu,label='NR-BWDEuler')  
   axs.plot(rkt,rku,label='Runge Kutta 4th')
   
   legend = leg = axs.legend()
   legend.get_frame().set_linewidth(1.2)
   plt.setp(legend.get_texts(), color='#555555')
   axs.set(xlabel = 'time' ,ylabel='y(t)')
   axs.set_title('Lorentz Attractor',y=1.05)
   plt.savefig('../Results/systemLorentz.pdf')   
   plt.show()

 
if __name__ == '__main__': 
   main()

