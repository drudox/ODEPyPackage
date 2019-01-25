# -*- coding : utf-8 -*-

import numpy as np
from ... Rhs.rhs import Rhs
from abc import abstractmethod
from .. ABSolver.sode import Sode
import scipy.optimize as opt
from scipy.optimize import fsolve, broyden1, newton_krylov

#

class BDF(Sode):
   def __init__(self, dydt : Rhs) :
      self.dydt   = dydt
      self.dt     = (dydt.tf-dydt.t0)/dydt.n
      self.f      = dydt.f
      #self._print = _print 
      #self.file   = filename
   
   def plot(self,*args):
       self.parameter= [*args]
       plt.style.use(['mystyle', 'mystyle-paper' ,'mystyle-nb'])
       ax = plt.subplot(111)
       ax.plot(self.time,self.u)
       ax.set_title (self.parameter[0])
       ax.set_xlabel(self.parameter[1])
       ax.set_ylabel(self.parameter[2])
       plt.show()
   
   @abstractmethod
   def solve(self):
      pass         

#---------------------------------------------------------------------------------------------------------

class MEBDF(BDF):
    '''
        Implement the solver Modified Extended BDF (Backward differences)
    '''
    
    solved = False
        
    def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            super().__init__(dydt)
    
    def stage1(self,variable,tp1,ti,ui,theta,c):
        yp1 = variable 
        f1  = -yp1 + ui + self.dt*(theta* self.f(tp1,yp1)+(1-theta)*self.f(ti,ui))
        return f1
    
    def stage2(self,variable,tp2,tp1,yp1,theta,c):
        yp2 = variable
        f2  = -yp2 + yp1 + self.dt*(theta* self.f(tp2,yp2)+(1+theta)* self.f(tp1,yp1))
        return f2
    
    def stage3(self,variable, tp2,tp1,ti,yp2,yp1,ui,theta,c):
        ynp1 = variable 
        f3   = -ynp1 + ui + self.dt*((c-0.5)* self.f(tp2,yp2)+(3/2 - 2*c -theta)* self.f(tp1,yp1) + theta*self.f(tp1,ynp1) + c*self.f(ti,ui) )
        return f3

    def solve(self,theta = 0.8 , c = 0):

        self.time , self.u = self.dydt.createArray()
        print ('running MEBDF method ...')
        for i in range(len(self.time)-1):
            
            up1         = fsolve(self.stage1, self.u[i], args=(self.time[i+1],self.time[i],self.u[i],theta,c))
            up2         = fsolve(self.stage2, up1      , args=(self.time[i+1]+self.dt,self.time[i+1],  up1  ,theta,c))    
            self.u[i+1] = fsolve(self.stage3, up1      , args=(self.time[i+1]+self.dt,self.time[i+1],self.time[i],up2,up1,self.u[i],theta,c))
        print('... Done!')
        return self.time,self.u
        ImpRK4.solved = True

        if ImpRK4.solved and self.dydt.solution[0] != 0:
            self.dydt.evalError(self.u, 'Imp RK-4th', prints = False)
               
        if self.file != None:
            super().write2file()

        if self.save:
              return self.time,self.u

    def plot(self):
            if ImpRK4.solved:
               super().plot('ODE Solution using Runge Kutta 2nd order','time [s]', 'y(t)')
            else:
               print("Unsolved problem, call `solve` method before")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

 

