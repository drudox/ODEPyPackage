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

func01 = lambda t,u : -.013*u[0] - 1000*u[0]*u[2]
func02 = lambda t,u : -2500*u[1]*u[2]
func03 = lambda t,u :-.013*u[0]-1000*u[0]*u[2] - 2500* u[1]*u[2]
#
#-------------------------------------------------------------
func0 = np.array([func01,func02,func03])
#anal0 = np.array([func1a,func2a])
def main():
   
   
   start_time = datetime.now()
      
   y0 = np.array([1.,1.,0.]) 
   
   
   problem1 = rhs.Rhs(func0, 0.0, 50 , y0, 6000)
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
   
   #beuler   = euler.Backward(problem1, 'bwe.dat')
   #bt,bu    = beuler.solve()
 
#   midpoint    = centdiff.MidPointSolver(problem1 , 'mp.dat')
#   cdt,cdu     = midpoint.solve()

   #sieuler_p4 = euler.Implicit(problem1, 'bwe.dat')
   #siet,sieu   = bwdeuler_p4.solve()
   #implicitmidpoint    = centdiff.ImplicitMidPoint(problem1 , 'Imp_mp.dat')
   #impt,impu           = implicitmidpoint.solve()



#   cn_p2    = rungekutta.ImplicitCrankNicholson(problem1, 'cn.dat')
#   cnt,cnu  = cn_p2.solve()
#    
#   rk_p1     =  implicitrungekutta.ImpRK2(problem1, 'rk.dat')
#   rkt,rku   =  rk_p1.solve()
  
   #rk_p     =  implicitrungekutta.ImpRK4(problem1, 'rk.dat')
   #rkt,rku  =  rk_p.solve()
   
   #bdf       =  MEBDF(problem1,'bdf.dat')
   #bdft,bdfu =  bdf.solve(0.8,0.5)
 
   #radau    =  RadauIIA(problem1,'Radsys3.dat')
   #rdt,rdu  =  radau.solve()
   
   radau5    =  RadauIIA5th(problem1,'Radsys3.dat')
   rd5t,rd5u  =  radau5.solve()


   end_time = datetime.now()
   print('Duration: {}'.format(end_time - start_time))
      
   
   fig,axs = qualityplot.Fonts()(1,1)    #mkplot.PlotBase('Solution of Differential Equations','p1.pdf')
   
 #  axs.plot(fet,feu,label='Exp Euler')
   #axs.plot(iet,ieu,label='Imp Euler')
   #axs.plot(siet,sieu,label='Sem.Imp Euler') #, linestyle=':')
   #axs.plot(bt,bu,label='NR-BWDEuler') #, linestyle=':')
   #axs.plot(iet,ieu,label='NR-BWDEuler') #, linestyle=':')
   #axs.plot(bet,beu,label='Imp Euler') #, linestyle=':')
   axs.plot(rd5t,rd5u,label='RadauIIA 5th')
   #axs.plot(cnt,cnu,label='Imp Crank Nicholson',linestyle='-')
   #axs.plot(rkt,rku,label='Imp RK4') #,linestyle=':')
  # axs.plot(rk3t,rk3u,label='Imp RK 3rd') #,linestyle=':')
 #  axs.plot(am4t,am4u,label='Adams Moulton 4th') #,linestyle=':')
   #axs.plot(bdft,bdfu,label='BDF')
   #axs.plot(impt,impu,label='Implicit Mid-Point Method',linestyle='-.')
   #axs.plot(x,y,label='Analitical', marker='o',markersize=7,linestyle='',alpha=0.7, color= 'C0' ) #, linestyle=':')
   #axs.plot(x,1.5*np.exp(-x) - 0.5* np.exp(-1000*x), marker='o',markersize=7,linestyle='',alpha=0.7 , color= 'C0') #, linestyle=':')
   #legend = leg = axs.legend()
   #legend.get_frame().set_linewidth(0.7)
   #plt.setp(legend.get_texts(), color='#555555')
   axs.set(xlabel = 'time' ,ylabel='y(t)')
   axs.set_title('Stiff Sys ODE',y=1.05)
    #plt.savefig('system02.pdf')   
   chartBox = axs.get_position()    
   legend = leg = axs.legend()
   legend.get_frame().set_linewidth(0.7)
   axs.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
   plt.setp(legend.get_texts(), color='#2E3436')
   plt.legend(bbox_to_anchor=(1.03, 0.75), loc=2, borderaxespad=0.)
  
   
   
   plt.show()
 
   
if __name__ == '__main__': 
   main()

