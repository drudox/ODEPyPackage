#!/usr/bin/env python3

'''   
  Driver programs for testing the solver 
  Sistem of first order differential equations 
  
  @author Marco Ghiani Glasgow 2018

  developing for python v.3.x 
'''


import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from ODE.Rhs import rhs
from ODE.qualityPlot import qualityplot
from ODE.Solvers.Euler import euler
from ODE.Solvers.RungeKutta import rungekutta
from ODE.Solvers.MultiStep import multistep
from ODE.Solvers.AdamsMethods import adamsmethods
from ODE.Solvers.CentDiff import centdiff
from ODE.Solvers.Euler import impliciteuler
from ODE.Solvers.RungeKutta import implicitrungekutta
from ODE.Solvers.AdamsMethods import adamsmoulton
from ODE.Solvers.RadauIIA.radauIIA import RadauIIA, RadauIIA5th
from ODE.Solvers.BDF.bdf import MEBDF



#------------------------------------------------------------
# Problem 1
#func01 = lambda t,u : -50*(u[0] - np.cos(t)) 
#anal01 = lambda x,y : np.exp(-5*(x-1)**2)

func01 = lambda t,u : -0.04* u[0]+ 10e4 * u[1] *(t-10e-2)*u[2] 
func02 = lambda t,u :  0.04* u[0]- 10e4 * u[1] *(t-10e-2)*u[2] - 3*10e7* u[1]**2
func03 = lambda t,u :  3 * 10e7 * u[1]**2

#
#-------------------------------------------------------------
## Lorentz attractor system

#func1a = lambda x,y : 2*np.exp(-x) - np.exp(-100*x)
#func2a = lambda x,y : 2*np.exp(-x) + np.exp(-100*x)

func0 = np.array([func01,func02,func03])
#anal0 = np.array([func1a,func2a])
def main():
   
   
   start_time = datetime.now()
      
   y0 = np.array([1.,0.,0]) 
   
   
   problem1 = rhs.Rhs(func0, 0.0, 10e11 , y0, 100000000)
   #x,y = problem1.analiticalSolution(save=True)   
   
 
##----------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------


#   fwdeuler_p1 = euler.Explicit(problem1,'fwe.dat')
#   fet,feu     = fwdeuler_p1.solve()
   
   #bwdeuler = impliciteuler.BWDEuler(problem1, 'bwe.dat')
   #iet,ieu  = bwdeuler.solve()


#   sieuler_p1  = euler.SemiImplicit(problem1, 'bwe.dat')
#   siet,sieu   = sieuler_p1.solve()
   
   beuler   = euler.Backward(problem1, 'bwe.dat')
   bt,bu    = beuler.solve()
 
#   midpoint    = centdiff.MidPointSolver(problem1 , 'mp.dat')
#   cdt,cdu     = midpoint.solve()

   #sieuler_p4 = euler.Implicit(problem1, 'bwe.dat')
   #siet,sieu   = bwdeuler_p4.solve()
#   implicitmidpoint    = centdiff.ImplicitMidPoint(problem1 , 'Imp_mp.dat')
#   impt,impu           = implicitmidpoint.solve()



#   cn_p2    = rungekutta.ImplicitCrankNicholson(problem1, 'cn.dat')
#   cnt,cnu  = cn_p2.solve()
#    
#   rk_p1     =  rungekutta.ImplicitRK2(problem1, 'rk.dat')
#   rkt,rku   =  rk_p1.solve()
  
   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
   
   
   fig,axs = qualityplot.Fonts()(1,1)    #mkplot.PlotBase('Solution of Differential Equations','p1.pdf')
   
 #  axs.plot(fet,feu,label='Exp Euler')
    #axs.plot(iet,ieu,label='Imp Euler')
   #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
   axs.plot(bt,bu,label='NR-BWDEuler') #, linestyle=':')
   #axs.plot(iet,ieu,label='NR-BWDEuler') #, linestyle=':')
   #axs.plot(bet,beu,label='Imp Euler') #, linestyle=':')
   #axs.plot(lft,lfu,label='Leap Frog')
   #axs.plot(cnt,cnu,label='Imp Crank Nicholson',linestyle='-')
   #axs.plot(rkt,rku,label='Imp RK2') #,linestyle=':')
   #axs.plot(am3t,am3u,label='Adams Moulton 3rd') #,linestyle=':')
 #  axs.plot(am4t,am4u,label='Adams Moulton 4th') #,linestyle=':')
   #axs.plot(cd2t,cd2u,label='Cent diff 2nd')
   #axs.plot(impt,impu,label='Implicit Mid-Point Method',linestyle='-.')
  # axs.plot(x,y,label='Analitical', marker='o',markersize=7,linestyle='',alpha=0.7, color= 'C0' ) #, linestyle=':')
   #axs.plot(x,1.5*np.exp(-x) - 0.5* np.exp(-1000*x), marker='o',markersize=7,linestyle='',alpha=0.7 , color= 'C0') #, linestyle=':')
   legend = leg = axs.legend()
   legend.get_frame().set_linewidth(0.7)
   plt.setp(legend.get_texts(), color='#555555')
   plt.xscale('log')
   axs.set(xlabel = 'time' ,ylabel='y(t)')
   axs.set_title('dy$_1$/dt=-10(t-1)y$_1$',y=1.05)
    #plt.savefig('system02.pdf')   
   plt.show()
 
   
if __name__ == '__main__': 
   main()

