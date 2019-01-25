#!/usr/bin/env python3

'''   
  Driver programs for testing the solver 
  Sistem of first order differential equations 
  
  @author Marco Ghiani Glasgow 2018

  developing for python v.3.x 
'''



import matplotlib.pyplot as plt
import numpy as np
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
func01 = lambda t,u : -1000*u[0] + 3000 - 2000* np.exp(-t) 
an01 = lambda x,y : 3 - 0.998*np.exp(-1000*x)-2.002*np.exp(-x)


func0 = np.array([func01])
anal01= np.array([an01])
def main():
   
   
   start_time = datetime.now()
      
   y0 = np.array([0.])
   
   problem1 = rhs.Rhs(func0, 0.0, 2.0 , y0, 30, anal01)
   
   x,y = problem1.analiticalSolution(save=True)   
 
#------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------
 
  
   bwdeuler    = euler.Implicit(problem1, 'bwe.dat')
   bet,beu     = bwdeuler.solve()


#   sieuler_p1  = euler.SemiImplicit(problem1, 'bwe.dat')
#   siet,sieu   = sieuler_p1.solve()
   
   beuler      = euler.Backward(problem1, 'bwe.dat')
   bt,bu       = beuler.solve()
  
   ieuler      = impliciteuler.BWDEuler(problem1, 'bwe.dat')
   iet,ieu     = ieuler.solve()
 
#   midpoint    = centdiff.MidPointSolver(problem1 , 'mp.dat')
#   cdt,cdu     = midpoint.solve()

   #sieuler_p4 = euler.Implicit(problem1, 'bwe.dat')
   #siet,sieu   = bwdeuler_p4.solve()
   implicitmidpoint    = centdiff.ImplicitMidPoint(problem1 , 'Imp_mp.dat')
   impt,impu           = implicitmidpoint.solve()

#   impmidpoint    = centdiff.ImpMidPoint(problem1 , 'Imp_mp.dat')
#   it,iu          = impmidpoint.solve()



   cn_p2    = rungekutta.ImplicitCrankNicholson(problem1, 'cn.dat')
   cnt,cnu  = cn_p2.solve()
#    
   rk_p1     =  implicitrungekutta.ImpRK2(problem1, 'rk.dat')
   rkt,rku  =  rk_p1.solve()
#
#   #fwdeuler_p1 = euler.Explicit(problem1, True, True, 'lorentz_fwe.dat')
#   #fet,feu = fwdeuler_p1.solve()
#   
   #leapfrog_p1 = multistep.LeapFrog(problem1, 'lorentz_lf.dat')
   #lft,lfu     = leapfrog_p1.solve()
#  
#   leapfrog_p2 = euler.Explicit(problem2, 'oscillator_lf.dat')
#   lf2,lf2     = leapfrog_p2.solve()


   rk       =  implicitrungekutta.ImpRK4(problem1, 'rk.dat')
   rk3t,rk3u  =  rk.solve()

   bdf       = MEBDF(problem1,'bdf.dat')
   bdft,bdfu = bdf.solve(1,0.5)

   radau     = RadauIIA(problem1,'radau.dat')
   rdt,rdu   = radau.solve()
   
   radau     = RadauIIA5th(problem1,'radau.dat')
   rd5t,rd5u   = radau.solve()
   
   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
   
   
   fig,axs = qualityplot.Fonts()(1,1)    #mkplot.PlotBase('Solution of Differential Equations','p1.pdf')
   
 #  axs.plot(fet,feu,label='Exp Euler')
   axs.plot(iet,ieu,label='Imp Euler')
   #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
   #axs.plot(bt,bu,label='NR-BWDEuler') #, linestyle=':')
 #  axs.plot(bet,beu,label='Imp Euler') #, linestyle=':')
   axs.plot(rdt,rdu,label='RadauIIA')
   axs.plot(rd5t,rd5u,label='RadauIIA 5th')
   #axs.plot(cnt,cnu,label='Implicit Crank Nicholson',linestyle='-')
   #axs.plot(rkt,rku,label='Runge Kutta 4th') #,linestyle=':')
   #axs.plot(am3t,am3u,label='Adams Moulton 3rd') #,linestyle=':')
   #axs.plot(am4t,am4u,label='Adams Moulton 4th') #,linestyle=':')
   axs.plot(bdft,bdfu,label='BDF')
 #  axs.plot(rkt,rku,label='Imp. Runge Kutta')
   #axs.plot(impt,impu,label='Implicit Mid-Point Method',linestyle='--')
   #axs.plot(iet,ieu,label='NR Implicit Euler') #,linestyle='-.')
   axs.plot(rk3t,rk3u,label='Imp. RK4') #,linestyle='-.')
   axs.plot(x,y,label='Analitical', marker='o',markersize='7',linestyle='',color='C0',alpha=0.7 ) #, linestyle=':')
   legend = leg = axs.legend()
   legend.get_frame().set_linewidth(0.7)
   plt.setp(legend.get_texts(), color='#555555')
   axs.set(xlabel = 'time' ,ylabel='y(t)')
   axs.set_title('dy$_1$/dt=-10(t-1)y$_1$',y=1.05)
    #plt.savefig('system02.pdf')   
   plt.show()
 
   
if __name__ == '__main__': 
   main()

