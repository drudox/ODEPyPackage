#!/usr/bin/env python3
'''   
   Multi-Step module 
   containers for classes used to solve System of DE 
   using implicit / Explicit multi-step methods 

'''

#from abc import ABCMeta , abstractmethod
from abc import abstractmethod
from .. ABSolver.sode import Sode
import numpy as np
from ... Rhs.rhs import Rhs
import matplotlib.pyplot as plt
import types
from .. RungeKutta import rungekutta

#---------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

class MultiStep(Sode):
   '''   
         Base class for all the MultiStep method 
         
   '''
   def __init__(self, dydt : Rhs) :
      self.dydt   = dydt
      self.dt     = (dydt.tf-dydt.t0)/dydt.n
      self.f      = dydt.f
      #self._print = _print 
      #self.file   = filename
   
   def plot(self,*args):
       self.parameter= [*args]
       plt.style.use(['mystyle', 'mystyle-nb'])
       ax = plt.subplot(111)
       ax.plot(self.time,self.u)
       ax.set_title (self.parameter[0])
       ax.set_xlabel(self.parameter[1])
       ax.set_ylabel(self.parameter[2])
       plt.show()
   
   @abstractmethod
   def solve(self):
      pass         


#-----------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#       Leap Frog (mid point) method
#---------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

class LeapFrog(MultiStep):
    '''
        perform the solution of a SODE using 2 step 
    '''
    solved = False
    um1 = 0.0
    startup = True

    def __init__(self, dydt : Rhs, filename: str=None, save: bool= True ):
        
        self.save = save
        self.file   = filename
        super().__init__(dydt) 
    
    def solve(self):
        ''' 
            Perform the leap-frog (mid point scheme)
        '''
        self.time , self.u = self.dydt.createArray()

        print("Running Leap-Frog (mid-point) ....")
        
        # start-up the method
        #self.u[1] = rungekutta.RK2.Heun.step(self.dydt.f, self.time[0], self.u[0], self.dt )
        k1 = self.dydt.f(self.time[0],self.u[0])
        k2 = self.dydt.f(self.time[0]+ self.dt , self.u[0]+ self.dt*k1 )

        self.u[1] = self.u[0] + self.dt/2.*(k1+k2)

        for i in range(1,len(self.time)-1):
            self.u[i+1] = self.u[i-1] + 2. * self.dt * self.dydt.f(self.time[i],self.u[i])
        
        print('... Done!')
        LeapFrog.solved = True 
        
        if self.save:
          return self.time,self.u

        if self.file != None:
           super().write2file() 
    
    def plot(self):
        if LeapFrog.solved:
           super().plot('ODE Solution using Leap Frog method','time [s]', 'y(t)')
        else:
           print("Unsolved problem, call `solve` method before")

#-----------------------------------------------------------------------------------------------------------

    @classmethod
    def step(cls ,func , t : np.float , u : np.float , dt ):
       
       def f(ti,ui):
            return  np.array([function(ti,ui) for function in func])     
       
  
       if cls.startup:
        
        
          k1 = f(t,u)
          k2 = f(t+dt, u+ dt*k1)
        
          unext   = u + dt/2.*(k1+k2)
          cls.startup = False
          cls.um1 = u.copy()    # u[i-1]
       else:
          unext = cls.um1 + dt *2.* f(t,u)
          cls.um1 = u
       return unext 
 






