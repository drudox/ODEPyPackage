#!/usr/bin/env python3
'''   
   Runge Kutta module 
   containers for classes used to solve System of DE 
   using implicit / Explicit RK methods 

'''

from abc import abstractmethod
from .. ABSolver.sode import Sode
from ... Rhs.rhs import Rhs
import numpy as np
import matplotlib.pyplot as plt
import types
from .. Euler import euler
#from Euler import *
import scipy.optimize as opt
from scipy.optimize import fsolve, broyden1, newton_krylov

#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

class RungeKutta(Sode):
   '''   
         Base class for all the Runge Kutta method 
         
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


#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------

class RK2(RungeKutta):
    '''
        solve a System Of Differential Equation (SODE) in time 
        using several Runge Kutta second order Method
    '''
    def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
        self.Heun(dydt,filename,save)
        self.ModifiedEuler(dydt,filename,save)
    
    @abstractmethod
    def solve(self):
        pass
#----------------------------------         Heun Method    -----------------------------------------------
        
    class Heun(RungeKutta):
        
        solved = False
        
        def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            super().__init__(dydt)

        def solve(self):
            '''
                perform the second order Heun scheme (RK)
            '''
            self.time , self.u = self.dydt.createArray()
        
            print("Running Runge-Kutta 2nd order (Heun's method) ....")
            for i in range(len(self.time)-1):
                k1 = self.dydt.f(self.time[i], self.u[i])
                k2 = self.dydt.f(self.time[i] + self.dt, self.u[i] + self.dt *k1 )
                self.u[i+1] = self.u[i] + self.dt/2.* (k1 + k2)      
            
            print('... Done!')
            RK2.Heun.solved = True
            
            if RK2.Heun.solved and self.dydt.solution[0] != 0:
                self.dydt.evalError(self.u, 'RK-2nd (Heun)', prints = False)
                 
            if self.file != None:
                super().write2file()


            if self.save:
                return self.time,self.u

            def plot(self):
                if RK2.Heun.solved:
                    super().plot('ODE Solution using Runge Kutta 2nd order','time [s]', 'y(t)')
                else:
                    print("Unsolved problem, call `solve` method before")
#
#---------------------------------------------------------------------------------------------------------
#
        @classmethod
        def step(cls,func , t : np.float , u : np.float , dt ):
        
            def f(ti,ui):
                return  np.array([function(ti,ui) for function in func])     
         
            k1 = f(t,u)
            k2 = f(t+dt , u + dt*k1)
            
            unext = u + dt/2.*(k1+k2)       
            return unext
        

#----------------------------------------------------------------------------------------------------------

#----------------------------------------  Modified Euler Method  -----------------------------------------

    class ModifiedEuler(RungeKutta):
       
       solved = False 
       
       def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            super().__init__(dydt)

       def solve(self):
            '''
                perform the second order Modified Euler (RK)

            '''
            self.time , self.u = self.dydt.createArray()
            print("Running Runge Kutta 2nd order (Modifier Euler's method) ....")
        
            for i in range(len(self.time)-1):
                k1 = self.dydt.f(self.time[i] ,  self.u[i])
                k2 = self.dydt.f(self.time[i] + self.dt/2. , self.u[i] + self.dt/2. *k1 )
                self.u[i+1] = self.u[i] + self.dt * k2      
            
            ('... Done!')
            RK2.ModifiedEuler.solved = True
            
            if RK2.ModifiedEuler.solved and self.dydt.solution[0] != 0:
                self.dydt.evalError(self.u, 'RK-2nd (Improved Euler)', prints = False)
            
            if self.file != None:
                super().write2file()

            if self.save:
                return self.time,self.u

            def plot(self):
                if RK2.ModifiedEuler.solved:
                    super().plot('ODE Solution using Runge Kutta 2nd order','time [s]', 'y(t)')
                else:
                    print("Unsolved problem, call `solve` method before")
#
#---------------------------------------------------------------------------------------------------------
#
       @classmethod
       def step(cls,func , t : np.float , u : np.float , dt ):
        
            def f(ti,ui):
                return  np.array([function(ti,ui) for function in func])     
         
            k1 = f(t,u)
            k2 = f(t + dt/2. , u + dt/2.*k1)
            
            unext = u + dt * k2       
            return unext
 

#----------------------------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------------------------
class CrankNicholson(RungeKutta):
   '''
      Cranck Nicholson system of differential equation solver
      Implicit Runge Kutta scheme
   '''
   solved = False

   def __init__(self, dydt : Rhs, filename: str=None, save :bool=True):
      
      self.save     = save
      self.file     = filename
      super().__init__(dydt)


   def solve(self):
      
      print("Running 2nd order Semi-Implicit (P-C) Crank-Nicholson method ....")
      self.time,self.u = self.dydt.createArray()
      for i in range(len(self.u)-1):
         '''
            Perform Predictor (Explicit Euler) - Corrector (2nd order Crank Niocholson)
         '''
         # Euler Predictor
         self.u[i+1] = self.u[i]+ self.dt* self.dydt.f(self.time[i],self.u[i])     
         # CrankNicholson corrector
         self.u[i+1] = self.u[i]+ self.dt/2 * ( self.dydt.f(self.time[i],self.u[i]) + self.dydt.f(self.time[i+1],self.u[i+1])) 
         
      print('... Done!')   
      CrankNicholson.solved=True
  
      if CrankNicholson.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Crank-Nicholson (Implicit-PC)', prints = False)
      
      if self.file != None:
          super().write2file()


      if self.save:
          return self.time,self.u

   def plot(self):
      if CrankNicholson.solved:
         super().plot('ODE Solution using Crank-Nicholson method','time [s]', 'y(t)')
      else:
         print("Unsolved problem, call `solve` method before")

#-----------------------------------------------------------------------------------------------------------
   
   @classmethod
   def step(cls,func , t : np.float , u : np.float , dt ):
       
       def f(ti,ui):
            return  np.array([function(ti,ui) for function in func])     
       
       # predictor
       up = euler.Explicit.step(func,t,u,dt) 
       # corrector 
       uc = u + dt/2 * (f(t,u) + f(t+dt,up) ) 
       unext = uc
       return unext

##----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
class ImplicitCrankNicholson(RungeKutta):
   '''
      Cranck Nicholson system of differential equation solver
      Implicit Runge Kutta scheme
   '''
   solved = False

   def __init__(self, dydt : Rhs, filename: str=None, save :bool=True):
      
      self.save     = save
      self.file     = filename
      super().__init__(dydt)


   def solve(self):
      
      toll = 10e-5
      print("Running 2nd order Implicit Crank-Nicholson method ....")
      self.time,self.u = self.dydt.createArray()
      for i in range(len(self.time)-1):
          
          func = lambda up1 : up1 - (self.u[i] + self.dt/2 *  (self.dydt.f(self.time[i],self.u[i])  +  self.dydt.f(self.time[i+1],up1)))
          
          self.u[i+1] = fsolve(func, self.u[i], xtol=10e-12)
         
      print('... Done!')   
      ImplicitCrankNicholson.solved=True
  
      if ImplicitCrankNicholson.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Crank-Nicholson (Implicit)', prints = False)
      
      if self.file != None:
          super().write2file()


      if self.save:
          return self.time,self.u

   def plot(self):
      if ImplicitCrankNicholson.solved:
         super().plot('ODE Solution using Implicit Crank-Nicholson method','time [s]', 'y(t)')
      else:
         print("Unsolved problem, call `solve` method before")

#-----------------------------------------------------------------------------------------------------------
 
class RK3(RungeKutta):
    
    solved = False

    def __init__(self, dydt : Rhs , filename: str = None , save : bool= True ):
        self.file = filename 
        self.save = save
        super().__init__(dydt)

    def solve(self):
        ''' 
            Solve the 3th order Runge Kutta Method
        '''
        print("Running Runge-Kutta 3rd order ....")
        self.time, self.u = self.dydt.createArray()

        for i in range(len(self.time)-1):
            k1 = self.dt * self.dydt.f(self.time[i],self.u[i])
            k2 = self.dt * self.dydt.f(self.time[i]+self.dt/2. , self.u[i]+ k1/2. )
            k3 = self.dt * self.dydt.f(self.time[i]+self.dt    , self.u[i]- k1 + 2*k2  )
            
            self.u[i+1] = self.u[i] + 1/6. * (k1 + 4*k2 + k3)

        print('... Done!')   
        RK3.solved = True
        
        if RK3.solved and self.dydt.solution[0] != 0:
           self.dydt.evalError(self.u, 'RK3 method', prints = False)
    
        if self.file != None:
            super().write2file()

        if self.save:
            return self.time , self.u
        
        def plot(self):
            if RK3.solved:
                super().plot('ODE Solution using Runge Kutta 3nd order','time [s]', 'y(t)')
            else:
                print("Unsolved problem, call `solve` method before")
##-----------------------------------------------------------------------------------------------------------

    @classmethod
    def step(cls,func , t : np.float , u : np.float , dt ):
       
       def f(ti,ui):
            return  np.array([function(ti,ui) for function in func])     
    
       k1 = dt* f(t,u)
       k2 = dt* f(t+dt/2.,u+k1/2.)
       k3 = dt* f(t+dt   , u - k1 +2*k2)
       
       unext = u + 1/6. * (k1 + 4*k2 + k3)
       return unext 
        


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
#         Runge Kutta  4th order (classic)
#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


class RK4(RungeKutta):
        
    solved = False
        
    def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            super().__init__(dydt)

    def solve(self):
        '''
        perform standard Fourth order RK scheme 
        '''
        self.time , self.u = self.dydt.createArray()
        print("Running Runge-Kutta 4th order") 
        for i in range(len(self.time)-1):
            k1 = self.dydt.f(self.time[i], self.u[i])
            k2 = self.dydt.f(self.time[i] + self.dt/2., self.u[i] + self.dt/2. *k1 )
            k3 = self.dydt.f(self.time[i] + self.dt/2., self.u[i] + self.dt/2. *k2 )
            k4 = self.dydt.f(self.time[i] + self.dt   , self.u[i] + self.dt    *k3 )       
            self.u[i+1] = self.u[i] + self.dt/6.* (k1 + 2*k2 + 2*k3 + k4)      
        print('... Done!')    
        RK4.solved = True
        
        if RK4.solved and self.dydt.solution[0] != 0:
           self.dydt.evalError(self.u, 'RK4 method', prints = False)
    
        if self.file != None:
            super().write2file()

        if self.save:
            return self.time,self.u

        def plot(self):
            if RK4.solved:
                super().plot('ODE Solution using Runge Kutta 4nd order','time [s]', 'y(t)')
            else:
                print("Unsolved problem, call `solve` method before")

#-----------------------------------------------------------------------------------------------------------

    @classmethod
    def step(cls , func ,t : np.float , u : np.float, dt ):
        
        def f(ti,ui):
            return  np.array([function(ti,ui) for function in func])     
        
        k1 = f(t,u)
        k2 = f(t + dt/2. , u + dt/2.*k1 )
        k3 = f(t + dt/2. , u + dt/2.*k2 )
        k4 = f(t + dt    , u + dt   *k3 )
        
        unext = u + dt/6. *( k1+ 2*k2 + 2*k3 + k4)       
        return unext


##-------------------------------------------------------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------------------------------------------------------
##                                         IMPLICIT RUNGE KUTTA METHODS  (STIFF SysODE SOLVER)
##-------------------------------------------------------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------------------------------------------------------

class ImplicitRK2(RungeKutta): ### ?????
        
    solved = False
        
    def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            super().__init__(dydt)

    def solve(self):
        '''
        perform standard Fourth order RK scheme 
        '''
        self.time , self.u = self.dydt.createArray()
        print("Running Implicit RK 2nd order") 
        for i in range(len(self.time)-1):
            
            func = lambda K1 : K1- self.f(self.time[i]+0.5*self.dt , self.u[i] + 0.5* K1) # ??
            
            k1 = newton_krylov(func , self.u[i] , f_tol = 10e-12 )                        # ??
            
            self.u[i+1] = self.u[i] + self.dt* k1                                         # ?? 
               
            
        print ('... Done!')

        ImplicitRK2.solved = True
        
        if ImplicitRK2.solved and self.dydt.solution[0] != 0:
            self.dydt.evalError(self.u, 'Implicit - Runge Kutta 2nd ord.', prints = False)
    
        if self.file != None:
            super().write2file()

        if self.save:
            return self.time,self.u

        def plot(self):
            if ImplicitRK2.solved:
                super().plot('ODE Solution using semi Implicit Runge Kutta Gill 4nd order','time [s]', 'y(t)')
            else:
                print("Unsolved problem, call `solve` method before")
 ##-------------------------------------------------------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------------
        
class RungeKuttaGill(RungeKutta):
        
        def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            super().__init__(dydt)

        def solve(self):
            '''
                perform the fourth order Runge Kutta Gill 
            '''
            self.time , self.u = self.dydt.createArray()
            a,b,c,d = (np.sqrt(2) - 1)/2. , (2 -np.sqrt(2))/2. , -np.sqrt(2)/2, (2+np.sqrt(2))/2. 
            
            print("Running Implicit Runge Kutta - Gill")       
            
            for i in range(len(self.time)-1):
                
                
                k1 = self.dt * self.f(self.time[i],self.u[i])
                k2 = self.dt * self.f(self.time[i] + 0.5*self.dt , self.u[i] + 0.5*k1) 
                k3 = self.dt * self.f(self.time[i] + 0.5*self.dt , self.u[i] + a*k1 + b*k2)
                k4 = self.dt * self.f(self.time[i] + self.dt     , self.u[i] + c*k2 + d*k3)           
                
                self.u[i+1] = self.u[i] + 1/6. *k1 + b/3. *k2 + d/3. *k3 + 1/6. *k4     
            
            print ('... Done!')

            RungeKuttaGill.solved = True
        
            if RungeKuttaGill.solved and self.dydt.solution[0] != 0:
                self.dydt.evalError(self.u, 'Semi Implicit - Runge Kutta Gill 4th ord.', prints = False)
    

            if self.file != None:
                super().write2file()

            if self.save:
                return self.time,self.u

            def plot(self):
                if RungeKuttaGill.solved:
                    super().plot('ODE Solution using semi Implicit Runge Kutta Gill 4nd order','time [s]', 'y(t)')
                else:
                    print("Unsolved problem, call `solve` method before")
    
##-------------------------------------------------------------------------------------------------------------------------------------------------------
