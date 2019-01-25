#!/usr/bin/env python3

'''
   Module to solve Ordinary differential equation 
   using first order Euler method (implicit and explicit)
   the pourpose of this module is to works with all 
   non-stiff differential problems 
   
   @author Marco Ghiani, Glasgow 2018
'''
try:
    from abc import ABCMeta,abstractmethod
    from .. ABSolver.sode import Sode
    from ... Rhs.rhs import Rhs
    import numpy as np
    import matplotlib.pyplot as plt
    import types 
except ImportError:
    print('Error occured importing modules!')
#--------------------------------------------------------------------------------------------

class CentrateDifference(Sode):
   ''' 
      Class base for the two Euler solvers (Explicit and Implicit)
   '''
   def __init__(self, dydt : Rhs) :
      self.dydt   = dydt
      self.f      = dydt.f
      self.dt     = (dydt.tf-dydt.t0)/dydt.n
   
   def plot(self,*args):
       '''  
         Plot the result of the sub classes integrations
         Implicit and Explicit method         
       ''' 
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
   
#---------------------------------------------------------------------------------------------
#        
#      
#
#
#---------------------------------------------------------------------------------------------

class MidPoint(CentrateDifference):
   '''
      Solve the ODE using first order centate difference
      mid point method
   '''
   
   solved = False  # True if the method 'solve' was called without error 

   def __init__(self, dydt: Rhs, filename : str=' ', save : bool=True):
      ''' 
         Initialize the midpoint method 
            - dydt (to the super class) : RHS problem dy/dt = f(t,y(t)) t(0)=y0 
            - save : if True returns the 2 vector time and u = du/dt 
      '''
      self.save   = save 
      self.file   = filename
    
      super().__init__(dydt) 

   def solve(self):
      self.time, self.u = self.dydt.createArray() 
      print('Running 2nd order Centrate Difference ....')
      for i in range(len(self.time)-1):
          
          self.du = self.dt* self.f(self.time[i],self.u[i]) 
          self.u[i+1] = self.u[i] + self.dt*self.f( self.time[i] + self.dt/2.  , self.u[i] + self.du/2. )
         
      print('... Done!')    
      MidPoint.solved = True 
      
      if MidPoint.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Cent Diff 2nd', prints = False)
      
      if self.file != None:
          super().write2file() 

      if self.save:
         return self.time,self.u
   
   def plot(self): 
     '''
         Check if the problem is solved (if the solve method was already called) 
         and if it is plot the numerical solution[0][0] of the IV differential problem
     '''
     if MidPoint.solved:
         super().plot('ODE Solution using Foward Euler', 'time [s]', 'y(t)')  
     else:    
         print("Unsolved problem, call `solve` method before")

#--------------------------------------------------------------------------------------------------

   @classmethod
   def step(cls,func , t : np.float , u : np.float , dt ):
       
       def f(ti,ui):
            return  np.array([function(ti,ui) for function in func])     
            
       du = dt*f(t,u) 
       unext = u + dt*f(t + dt/2. , u + du/2. )
       return unext


#-----------------------------------------------------------------------------------------------------
class   CentDiff2nd(CentrateDifference):
   '''
      Solve the ODE using first order centate difference
      mid point method
   '''
   
   solved = False  # True if the method 'solve' was called without error 

   def __init__(self, dydt: Rhs, filename : str=' ', save : bool=True):
      ''' 
         Initialize the midpoint method 
            - dydt (to the super class) : RHS problem dy/dt = f(t,y(t)) t(0)=y0 
            - save : if True returns the 2 vector time and u = du/dt 
      '''
      self.save   = save 
      self.file   = filename
    
      super().__init__(dydt) 

   def solve(self):
      self.time, self.u = self.dydt.createArray() 
      print('Running 2nd order Centrate Difference ....')
      for i in range(len(self.time)-1):
          
          self.u[i+1] = self.u[i] + self.dt*self.f( self.time[i] + self.dt/2.  , self.u[i] + self.dt/2.* self.f(self.time[i],self.u[i]) )
         
      print('... Done!')    
      MidPoint.solved = True 
      
      if MidPoint.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Cent Diff 2nd', prints = False)
      
      if self.file != None:
          super().write2file() 

      if self.save:
         return self.time,self.u
   
   def plot(self): 
     '''
         Check if the problem is solved (if the solve method was already called) 
         and if it is plot the numerical solution[0][0] of the IV differential problem
     '''
     if MidPoint.solved:
         super().plot('ODE Solution using Foward Euler', 'time [s]', 'y(t)')  
     else:    
         print("Unsolved problem, call `solve` method before")

#--------------------------------------------------------------------------------------------------

   @classmethod
   def step(cls,func , t : np.float , u : np.float , dt ):
       
       def f(ti,ui):
            return  np.array([function(ti,ui) for function in func])     
            
       #du = dt*f(t,u) 
       unext = u + dt*f(t + dt/2. , u + dt/2.*f(t,u) )
       return unext

#---------------------------------------------------------------------------------------------------
class MidPointSolver(CentrateDifference):
   '''
      Solve the ODE using first order centate difference
      mid point method
   '''
   
   solved = False  # True if the method 'solve' was called without error 
   startup = True 
   up = np.array([0])
   def __init__(self, dydt: Rhs, filename : str=' ', save : bool=True):
      ''' 
         Initialize the midpoint method 
            - dydt (to the super class) : RHS problem dy/dt = f(t,y(t)) t(0)=y0 
            - save : if True returns the 2 vector time and u = du/dt 
      '''
      self.save   = save 
      self.file   = filename
    
      super().__init__(dydt) 

   def solve(self):
      self.time, self.u = self.dydt.createArray() 
      print('Running 2nd order Mid Point Method ....')
        
      self.u[1] = self.u[0] + self.dt* self.f(self.time[0],self.u[0])

      for i in range(1,len(self.time)-1):
          self.u[i+1] = self.u[i-1] + 2*self.dt* self.f(self.time[i], self.u[i]) 
      
      print('... Done!')  
      MidPointSolver.solved = True 
      
      if MidPoint.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Cent Diff 2nd', prints = False)
      
      if self.file != None:
          super().write2file() 

      if self.save:
         return self.time,self.u
   
   def plot(self): 
     '''
         Check if the problem is solved (if the solve method was already called) 
         and if it is plot the numerical solution[0][0] of the IV differential problem
     '''
     if MidPointSolver.solved:
         super().plot('ODE Solution using Foward Euler', 'time [s]', 'y(t)')  
     else:    
         print("Unsolved problem, call `solve` method before")

#-------------------------------------------------------------------------------------------------------------------------------------------     
            

   @classmethod
   def step(cls,func , t : np.float , u : np.array , dt ):
       
        def f(ti,ui):
            return  np.array([function(ti,ui) for function in func])     
        
        if cls.startup:
            cls.up = u + dt* f(t,u)
            cls.startup = False
        unext = cls.up + 2*dt*f(t,u)
        cls.up = u.copy()
        return unext

#--------------------------------------------------------------------------------------------------------------------------------------------
#    Implicit midpoint method 
#---------------------------------------------------------------------------------------------------------------------------------------------


class ImplicitMidPoint(CentrateDifference):
   '''
      Solve the ODE using second order implicit mid point 
      Implicit mid-point method
   '''
   
   solved = False  # True if the method 'solve' was called without error 

   def __init__(self, dydt: Rhs, filename : str=' ', save : bool=True):
      ''' 
         Initialize the midpoint method 
            - dydt (to the super class) : RHS problem dy/dt = f(t,y(t)) t(0)=y0 
            - save : if True returns the 2 vector time and u = du/dt 
      '''
      self.save   = save 
      self.file   = filename
    
      super().__init__(dydt) 

   def solve(self):

      toll = 10e-6  
      self.time, self.u = self.dydt.createArray() 
      
      print('Running Implicit Mid-Point method....')
      
      for i in range(len(self.time)-1):
          
          # running implicit euler half step 

          u0 = self.u[i].T

          u1_2 = u0 - np.linalg.inv(np.eye(len(self.dydt.u0)) - self.dt/2. * self.dydt.df(self.time[i+1]- self.dt/2. , u0)).dot( u0 - self.dt/2. *self.f(self.time[i+1]-self.dt/2. ,u0).T  - self.u[i].T )

          error = np.array([1.0])
          iters = 0 
          while True:
              try:
                  u0= u1_2.T
                  u1_2 = u0 - np.linalg.inv(np.eye(len(self.dydt.u0)) - self.dt/2. * self.dydt.df(self.time[i+1]- self.dt/2. , u0)).dot( u0 - self.dt/2. *self.f(self.time[i+1]-self.dt/2. ,u0).T  - self.u[i].T )
                  
                  iters+=1 
                  error = np.abs(u1_2 - u0)
                  
                  if np.sum(np.abs(error)) <= toll :
                      break
                  if iters >= 1000: 
                      raise ValueError('too many iteractions')
              except ValueError as ex:
                  print('Implicit Mid Point, ERROR : %s' %ex.args)

          uHalf = u1_2.T              
              
          self.u[i+1] = uHalf + self.dt/2. * self.f(self.time[i]+ self.dt/2. , uHalf)
        
      print('... Done!') 
      ImplicitMidPoint.solved = True 
      
      if ImplicitMidPoint.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Implicit MidPoint', prints = False)
      
      if self.file != None:
          super().write2file() 

      if self.save:
         return self.time,self.u
   
   def plot(self): 
     '''
         Check if the problem is solved (if the solve method was already called) 
         and if it is plot the numerical solution[0][0] of the IV differential problem
     '''
     if ImplicitMidPoint.solved:
         super().plot('ODE Solution using Foward Euler', 'time [s]', 'y(t)')  
     else:    
         print("Unsolved problem, call `solve` method before")

#---------------------------------------------------------------------------------------------------------------------------------------------


class ImpMidPoint(CentrateDifference):
   '''
      Solve the ODE using second order implicit mid point 
      Implicit mid-point method
   '''
   
   solved = False  # True if the method 'solve' was called without error 

   def __init__(self, dydt: Rhs, filename : str=' ', save : bool=True):
      ''' 
         Initialize the midpoint method 
            - dydt (to the super class) : RHS problem dy/dt = f(t,y(t)) t(0)=y0 
            - save : if True returns the 2 vector time and u = du/dt 
      '''
      self.save   = save 
      self.file   = filename
    
      super().__init__(dydt) 

   def solve(self):

      toll = 10e-6  
      self.time, self.u = self.dydt.createArray() 
      
      print('Running Implicit Mid-Point method....')
      
      for i in range(len(self.time)-1):
          
          # running implicit euler half step 

          u0 = self.u[i].T

          u1_2 = u0 - np.linalg.inv(np.eye(len(self.dydt.u0)) - self.dt/2. * self.dydt.df(self.time[i+1]- self.dt/2. , u0)).dot( u0 - self.dt/2. *self.f(self.time[i+1]-self.dt/2. ,u0).T  - self.u[i].T )

          error = np.array([1.0])
          iters = 0 
          while True:
              try:
                  u0= u1_2.T
                  u1_2 = u0 - np.linalg.inv(np.eye(len(self.dydt.u0)) - self.dt/2. * self.dydt.df(self.time[i+1]- self.dt/2. , u0)).dot( u0 - self.dt/2. *self.f(self.time[i+1]-self.dt/2. ,u0).T  - self.u[i].T )
                  
                  iters+=1 
                  error = np.abs(u1_2 - u0)
                  
                  if np.sum(np.abs(error)) <= toll :
                      break
                  if iters >= 1000: 
                      raise ValueError('too many iteractions')
              except ValueError as ex:
                  print('NR Mid Point, ERROR : %s' %ex.args)

          u0 = u1_2.T

          u1 = u0 - np.linalg.inv(np.eye(len(self.dydt.u0)) - self.dt/2. * self.dydt.df(self.time[i+1] , u0)).dot( u0 - self.dt/2. *self.f(self.time[i+1] ,u0).T  - self.u[i].T )
          error = 1.0
          iters = 0 
       
          while True:
              try:
                  u0= u1.T
          
                  u1 = u0 - np.linalg.inv(np.eye(len(self.dydt.u0)) - self.dt/2. * self.dydt.df(self.time[i+1] , u0)).dot( u0 - self.dt/2. *self.f(self.time[i+1] ,u0).T  - self.u[i].T )
          
                  iters +=1 
                  error = np.abs(u1 - u0)
                  
                  if np.sum(np.abs(error)) <= toll :
                      break
                  if iters >= 1000: 
                      raise ValueError('too many iteractions')
              except ValueError as ex:
                  print('NR Mid Point, ERROR : %s' %ex.args)

          self.u[i+1] = u1.T #uHalf + self.dt/2. * self.f(self.time[i]+ self.dt/2. , uHalf)
        
      print('... Done!') 
      ImplicitMidPoint.solved = True 
      
      if ImplicitMidPoint.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Implicit MidPoint', prints = False)
      
      if self.file != None:
          super().write2file() 

      if self.save:
         return self.time,self.u
   
   def plot(self): 
     '''
         Check if the problem is solved (if the solve method was already called) 
         and if it is plot the numerical solution[0][0] of the IV differential problem
     '''
     if ImplicitMidPoint.solved:
         super().plot('ODE Solution using Foward Euler', 'time [s]', 'y(t)')  
     else:    
         print("Unsolved problem, call `solve` method before")


##------------------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------------------------------
 





