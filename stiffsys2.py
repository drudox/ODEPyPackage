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

func01 = lambda t,u : -500.5 * u[0] + 499.5 *u[1]
func02 = lambda t,u : 499.5*u[0] - 500.5* u[1]

#
#-------------------------------------------------------------
## Lorentz attractor system

func1a = lambda x,y : 1.5*np.exp(-x) + 0.5* np.exp(-1000*x)
func2a = lambda x,y : 1.5*np.exp(-x) - 0.5* np.exp(-1000*x)

func0 = np.array([func01,func02])
anal0 = np.array([func1a,func2a])
def main():
   
   
   start_time = datetime.now()
      
   #y0_sode_p1 = np.array([1.,0.,0.])
   y0 = np.array([2.,1.]) 
   
   
 #  problem1 = rhs.Rhs(func01, 0.0, 2 , y0 , 30)
   


   #y0_sode_p1 = np.array([1.,0.,0.])
   #y0_sode_p2 = np.array([1.,0.]) 
   
   #y0 = np.array([0.15])
   
   problem1 = rhs.Rhs(func0, 0.0, 0.01 , y0, 20, anal0)
   x,y = problem1.analiticalSolution(save=True)   
   
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------


#   fwdeuler_p1 = euler.Explicit(problem1,'fwe.dat')
#   fet,feu     = fwdeuler_p1.solve()
   
   bwdeuler    = impliciteuler.BWDEuler(problem1, 'bwe.dat')
   bet,beu     = bwdeuler.solve()


#   sieuler_p1  = euler.SemiImplicit(problem1, 'bwe.dat')
#   siet,sieu   = sieuler_p1.solve()
   
   beuler      = euler.Backward(problem1, 'bwe.dat')
   bt,bu       = beuler.solve()
 
#   midpoint    = centdiff.MidPointSolver(problem1 , 'mp.dat')
#   cdt,cdu     = midpoint.solve()

   #sieuler_p4 = euler.Implicit(problem1, 'bwe.dat')
   #siet,sieu   = bwdeuler_p4.solve()
   implicitmidpoint    = centdiff.ImplicitMidPoint(problem1 , 'Imp_mp.dat')
   impt,impu           = implicitmidpoint.solve()



   impcn      = rungekutta.ImplicitCrankNicholson(problem1, 'cn.dat')
   cnt,cnu  =  impcn.solve()
#    
#   rk_p1     =  implicitrungekutta.ImpRK2(problem1, 'rk.dat')
#   rkt,rku   =  rk_p1.solve()
   
   rk_p1       =  implicitrungekutta.ImpRK4(problem1, 'rk.dat')
   rkt,rku   =  rk_p1.solve()
   
   radau     = RadauIIA(problem1,'radIIA.dat')
   rdt,rdu   = radau.solve()
  
   radau5     = RadauIIA5th(problem1,'radIIA.dat')
   rd5t,rd5u  = radau5.solve()
   
   bdf       = MEBDF(problem1,'bdf.dat')
   bdft,bdfu = bdf.solve(1,0.5)

   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
   
   
   fig,axs = qualityplot.Fonts()(1,1)    #mkplot.PlotBase('Solution of Differential Equations','p1.pdf')
   
 #  axs.plot(fet,feu,label='Exp Euler')
    #axs.plot(iet,ieu,label='Imp Euler')
   #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
   axs.plot(bt,bu,label='NR-BWDEuler',color='C1') #, linestyle=':')
   #axs.plot(bet,beu,label='Imp Euler') #, linestyle=':')
   axs.plot(rdt,rdu,label='Radau IIA 3rd',color='C3')
   axs.plot(rd5t,rd5u,label='Radau IIA 5th',color='C2')
   axs.plot(bdft,bdfu,label='BDF',color='C4')#,linestyle='--')
   axs.plot(rkt,rku,label='Implicit RK4',color='C0') #,linestyle=':')
   #axs.plot(rk3t,rk3u,label='Imp RK') #,linestyle=':')
   #axs.plot(am4t,am4u,label='Adams Moulton 4th') #,linestyle=':')
   #axs.plot(cd2t,cd2u,label='Cent diff 2nd')
  # axs.plot(impt,impu,label='Implicit Mid-Point Method')#,linestyle='-.')
   axs.plot(x,1.5*np.exp(-x) + 0.5* np.exp(-1000*x) ,label='Analitical', marker='o',markersize=7,linestyle='',alpha=0.7, color= 'C0' ) #, linestyle=':')
   axs.plot(x,1.5*np.exp(-x) - 0.5* np.exp(-1000*x), marker='o',markersize=7,linestyle='',alpha=0.7 , color= 'C0') #, linestyle=':')
   legend = leg = axs.legend()
   legend.get_frame().set_linewidth(0.7)
   plt.setp(legend.get_texts(), color='#555555')
   axs.set(xlabel = 'time' ,ylabel='y(t)')
   axs.set_title('dy$_1$/dt=-10(t-1)y$_1$',y=1.05)
    #plt.savefig('system02.pdf')   
   chartBox = axs.get_position()    
   legend = leg = axs.legend()
   legend.get_frame().set_linewidth(0.7)
   axs.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
   plt.setp(legend.get_texts(), color='#2E3436')
   plt.legend(bbox_to_anchor=(1.03, 0.85), loc=2, borderaxespad=0.)
   
   plt.show()
   
if __name__ == '__main__': 
   main()

