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

func01 = lambda t,u : -10*u[0]+ 6*u[1] 
func02 = lambda t,u : 13.5*u[0] -10*u[1]

func0 = np.array([func01,func02])

def main():
   
   
   start_time = datetime.now()
      
   y0 = np.array([2./3. * np.exp(2), 1 ]) 
   
  
   problem1 = rhs.Rhs(func0, 0.0, 3.0 , y0, 30 )
   #x,y = problem1.analiticalSolution(save=True)   

##------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------
   
   fwdeuler_p1 = euler.Explicit(problem1,'fwe.dat')
   fet,feu     = fwdeuler_p1.solve()
   
   bwdeuler_p1 = euler.Implicit(problem1, 'bwe.dat')
   bet,beu     = bwdeuler_p1.solve()


#   sieuler_p1  = euler.SemiImplicit(problem1, 'bwe.dat')
#   siet,sieu   = sieuler_p1.solve()
   
   beuler_p1   = euler.Backward(problem1, 'bwe.dat')
   bt,bu       = beuler_p1.solve()
 
#   midpoint    = centdiff.MidPointSolver(problem1 , 'mp.dat')
#   cdt,cdu     = midpoint.solve()

   #sieuler_p4 = euler.Implicit(problem1, 'bwe.dat')
   #siet,sieu   = bwdeuler_p4.solve()
   implicitmidpoint    = centdiff.ImplicitMidPoint(problem1 , 'Imp_mp.dat')
   impt,impu           = implicitmidpoint.solve()
 
   am5_p1 = adamsmoulton.AdamsMoulton5th(problem1, 'am5_1.dat')
   #ab2_p  = ab2_p1._2nd(problem1,'ab2_1.dat')
   am5t,am5u = am5_p1.solve()
 


   impcn      = implicitrungekutta.ImpRK4(problem1, 'cn.dat')
   rkt,rku    = impcn.solve()
   
   radau      = RadauIIA(problem1, 'radIIA.dat')
   radt,radu  = radau.solve()
   
   radau5       = RadauIIA5th(problem1, 'radIIA.dat')
   rad5t,rad5u  = radau5.solve()



 #  rk_p1     =  rungekutta.RK4(problem1, 'rk.dat')
 #  rkt,rku  =  rk_p1.solve()
#
   #rk3        =  implicitrungekutta.ImpRK4(problem1, 'rk.dat')
   #rk3t,rk3u  =  rk3.solve()
#   #fwdeuler_p1 = euler.Explicit(problem1, True, True, 'lorentz_fwe.dat')
#   #fet,feu = fwdeuler_p1.solve()
#   
   #leapfrog_p1 = multistep.LeapFrog(problem1, 'lorentz_lf.dat')
   #lft,lfu     = leapfrog_p1.solve()
#  
#   leapfrog_p2 = euler.Explicit(problem2, 'oscillator_lf.dat')
#   lf2,lf2     = leapfrog_p2.solve()





   #fwdeuler_p1.plot()
   
   #ab2_p1 = euler.Implicit(problem3,True,True,filename='lorentz_bwe.dat')
   #bet,beu = bwdeuler_p1.solve()
   #bwdeuler_p1.plot()
   
   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
   
   
   fig,axs = qualityplot.Fonts()(1,1)    #mkplot.PlotBase('Solution of Differential Equations','p1.pdf')
   
 #  axs.plot(fet,feu,label='Exp Euler')
    #axs.plot(iet,ieu,label='Imp Euler')
   #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
   axs.plot(bt,bu,label='NR-BWDEuler') #, linestyle=':')
 #  axs.plot(bet,beu,label='Imp Euler') #, linestyle=':')
   axs.plot(radt,radu,label='RadauIIA 3th',linewidth=3)
   axs.plot(rad5t,rad5u,label='RadauIIA 5th',linewidth=3)
   axs.plot(rkt,rku,label='Imp RK4',linewidth=2)
   #axs.plot(am5t,am5u,label='Adams Moulton') #,linestyle=':')
   #axs.plot(am3t,am3u,label='Adams Moulton 3rd') #,linestyle=':')
   #axs.plot(am4t,am4u,label='Adams Moulton 4th') #,linestyle=':')
   #axs.plot(cd2t,cd2u,label='Cent diff 2nd')
   #axs.plot(impt,impu,label='Implicit Mid-Point Method',linestyle='-.')
   #axs.plot(x,y,label='Analitical', marker='o',linestyle='',color='C0',alpha=0.7 ) #, linestyle=':')
   legend = leg = axs.legend()
   legend.get_frame().set_linewidth(0.7)
   plt.setp(legend.get_texts(), color='#555555')
   axs.set(xlabel = 'time' ,ylabel='y(t)')
   axs.set_title('dy$_1$/dt=-10(t-1)y$_1$',y=1.05)
    #plt.savefig('system02.pdf')   
   plt.show()
 
   
if __name__ == '__main__': 
   main()

