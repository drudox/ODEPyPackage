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
    #from Solvers.ABSolver.sode import Sode
    from .. ABSolver.sode import Sode
    from ... Rhs.rhs import Rhs
    import numpy as np
    import matplotlib.pyplot as plt
    import types 
except ImportError:
    print('Error occured importing modules!')
#--------------------------------------------------------------------------------------------

class Euler(Sode):
   ''' 
      Class base for the two Euler solvers (Explicit and Implicit)
   '''
   def __init__(self, dydt : Rhs) :
      self.dydt   = dydt
      self.dt     = (dydt.tf-dydt.t0)/dydt.n
      self.f      = dydt.f    
   def plot(self,*args):
       '''  
         Plot the results of the sub classes integrations
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

class Explicit(Euler):
   '''
      Solve the ODE using first order Foward difference
      Explicit first order Method (Euler)  
   '''
   
   solved = False  # True if the method 'solve' was called without error 

   def __init__(self, dydt: Rhs, filename : str=' ', save : bool=True):

      ''' 
         Initialize the Foward Euler solver 
            - dydt (to the super class) : RHS problem dy/dt = f(t,y(t)) t(0)=y0 
            - save : if True returns the 2 vector time and u = du/dt 
      '''

      self.save   = save 
      self.file   = filename

      super().__init__(dydt) 

   def solve(self):
      self.time, self.u = self.dydt.createArray() 
      print('Running Explicit Euler ....')
      for i in range(len(self.time)-1):
              
          self.u[i+1] = self.u[i] + self.dt*self.dydt.f(self.time[i],self.u[i])
         
      print('... Done!')    
      Explicit.solved = True 
      
      if Explicit.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Explicit Euler', prints=False)

      if self.file != None:
          super().write2file() 

      if self.save:
         return self.time,self.u
    


   def plot(self): 
     '''
         Check if the problem is solved (if the solve method was already called) 
         and if it is plot the numerical solution[0] of the IV differential problem
     '''
     if Explicit.solved:
         super().plot('ODE Solution using Foward Euler', 'time [s]', 'y(t)')  
     else:    
         print("Unsolved problem, call `solve` method before")

#--------------------------------------------------------------------------------------------------

   @classmethod
   def step(cls,func , t : np.float , u : np.float , dt ):
       
       def f(ti,ui):
            return  np.array([function(ti,ui) for function in func])     
        
       unext = u + dt*f(t,u)
       return unext



  
#---------------------------------------------------------------------------------------------------
#
#
#---------------------------------------------------------------------------------------------------

class Implicit(Euler):
   '''
      Solve the ODE using first order Backward difference
      Implicit Method (Euler)  
   '''
   solved = False
   def __init__(self, dydt: Rhs, filename : str=' ', save :bool=True):
      ''' 
         Initialize the Implicit Euler solver 
            - dydt : RHS problem dy/dt = f(t,y(t)) t(0)=y0 
            - save : if True returns the 2 vector time and u = du/dt 
      '''
      self.save   = save
      self.file   = filename
      super().__init__(dydt)

   def solve(self):
      '''
            Solve the non linear implicit euler equations 
            using Newton Raphson method for the NLE 
            this method have a very large region of convergence 
            
            @Marco Ghiani PhD Glasgow 
      '''

      print('Running Implicit Euler ....')
      toll = 10e-6
      self.time, self.u = self.dydt.createArray()
      for i in range (len(self.time)-1):
         
         u_old = self.u[i] + self.dydt.f(self.time[i+1],self.u[i])
         error = np.array([1.0])
         
         iters = 0
         while True:
            try:            
                u_new = u_old - ((u_old - (self.u[i]+ self.dt*self.dydt.f(self.time[i+1], u_old)))/
                                  (1 - self.dt*self.dydt.Df(self.time[i+1],u_old) ) )
                error = np.abs(u_new - u_old)
                u_old = u_new
                iters+=1
                if np.sum(np.abs(error)) <= toll :
                 break
                if iters >= 1000:
                     raise ValueError('Too many iteraction')
            except ValueError as ex:
                print('Error in Implicit Euler Solver %s',ex)
                return None,None
            #else:
         self.u[i+1] = u_new

      print('... Done!')  
      Implicit.solved = True 
      
      if Implicit.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Implicit Euler', False)

 
      if self.file != None:
          super().write2file()

        #  with open(self.file,'w') as f:
        #      for i in range(len(self.time)):
        #          f.write('%.4f ' %self.time[i])
        #          for j in range(len(self.u[0])): 
        #              f.write('  %.4f  ' %self.u[i,j])        #, self.u[i,1], self.u[i,2]))
        #          f.write('\n')  
 
 
      if self.save:
        return self.time,self.u


   def plot(self):
     '''
         Check if the problem is solved (if the solve method was already called) 
         and if it is plot the numerical solution[0] of the IV differential problem
     '''
     if Implicit.solved:
         super().plot('ODE solution[0] using Backward Euler', 'time [s]', 'y(t)')    
     else:
         print("Unsolved problem, call `solve` method before")

#----------------------------------------------------------------------------------------------------------
# u_new = u_old - ((u_old - (self.u[i]+ self.dt*self.dydt.f(self.time[i+1], u_old)))/
#                              (1 - self.dt*self.dydt.Df(self.time[i+1],u_old) ) )
#            error = np.abs(u_new[:] - u_old[:])
#            u_old = u_new
#            if error.all() <= toll :
#               break
         
#            self.u[i+1] = u_new



   @classmethod
   def step(cls ,func , t : np.float , u : np.array , dt ):
       '''
            Implicit Euler Stepper, perform one step a time 
       '''
       def f(t,u):
            return  np.array([function(t,u) for function in func])     
       def Df(t,u):
            eps = 10e-12
            return ((f(t,u)+eps) - f(t,u))/eps
        
       uold  = Explicit.step(func,t+dt,u,dt) 
       
       toll  = 1e-6
       error = np.array([1.0])
       
       iters=0
       while True:  # Newton Rapson solving the non linearity
         try:  
           unew = uold - ( (uold - (u + dt*f(t+dt,uold) ) ) / (1 - dt*Df(t+dt,uold) ) )
           
           error = np.abs(unew - uold)
           uold  = unew 
           iters +=1 
           if np.sum(np.abs(error)) <= toll :
                break
           if iters >= 1000:
                raise ValueError('Too many Iteractions')
         except ValueError as ex:
                print('Error occurred in Implicit Euler Stepper %s', ex.args)
                return None
       return unew
##########################################################################################################
#                       Semi Implicit (predictor corrector) Euler solver
##########################################################################################################
class SemiImplicit(Euler):
   '''
      Solve the ODE using first order predictor-corrector 
      Semi-Implicit Method (Euler)  
   '''
   solved = False
   def __init__(self, dydt: Rhs, filename : str=' ', save :bool=True):
      ''' 
         Initialize the Semi Implicit Euler solver 
            - dydt : RHS problem dy/dt = f(t,y(t)) t(0)=y0 
            - save : if True returns the 2 vector time and u = du/dt 
      '''
      self.save   = save
      self.file   = filename
      super().__init__(dydt)

   def solve(self):
      '''
            Solve the non linear implicit euler equations 
            using Newton Raphson method for the NLE 
            this method have a very large region of convergence 
            
            @Marco Ghiani PhD Glasgow 
      '''
      print('Running Semi Implicit Euler ....')
      toll = 10e-8
      self.time, self.u = self.dydt.createArray()
      for i in range (len(self.time)-1):
            up = self.u[i] + self.dt * self.dydt.f(self.time[i+1], self.u[i])
            
            iter=0
            while True:
                uc = self.u[i] + self.dt * self.dydt.f(self.time[i+1], up)
                up = uc
                iter +=1
                if(iter >= 0):
                    break
            self.u[i+1] = uc
      
      print('... Done!')
      SemiImplicit.solved = True 
      
      if SemiImplicit.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Semi-Imp. Euler', prints = False)


      if self.file != None:
          super().write2file()

      if self.save:
        return self.time,self.u


   def plot(self):
     '''
         Check if the problem is solved (if the solve method was already called) 
         and if it is plot the numerical solution[0] of the IV differential problem
     '''
     if SemiImplicit.solved:
         super().plot('ODE solution[0] using Backward Euler', 'time [s]', 'y(t)')    
     else:
         print("Unsolved problem, call `solve` method before")
#--------------------------------------------------------------------------------------------------------
#    stepper 
   
   @classmethod 
   def step(cls , func , t :'np.float' , u : 'np.array', dt):
        
        def f(t,u):
            return np.array([function(t,u) for function in func])
        
        # explicit predictor
        up = Explicit.step(func,t+dt,u,dt)
        it=0
        while True:
            uc = u + dt*f(t+dt,up) 
            up = uc
            it +=1
            if (it >= 1):
                break
        return uc

#########################################################################################################
#                                   pure Implicit Backward Euler 
#########################################################################################################

class Backward(Euler):
    ''' 
        Backward Euler , Implicit scheme .
        Newton-Raphson used with the purpose to solve the non linear Equation
    '''
    solved = False
    def __init__(self, dydt: Rhs, filename : str=' ', save :bool=True):
       ''' 
            Initialize pure Implicit Euler solver 
            - dydt : RHS problem dy/dt = f(t,y(t)) t(0)=y0 
            - save : if True returns the 2 numpy.array (time and u = du/dt)
            - filename : if not ' ', writes the results in a file named filename 
       '''
       self.save   = save
       self.file   = filename
       super().__init__(dydt)

    def solve(self):
        '''
              Solve the non linear Implicit Euler solver 
              using Newton Raphson method for solve the NLE (non linear equation) each step
              this method have a very large region of convergence 
              
              @Marco Ghiani PhD Glasgow 
        '''
        toll = 10e-6
        self.time, self.u = self.dydt.createArray()   
        
        print('Running Newton-Rapson Backward Euler ....')
        
        for i in range(len(self.time)-1):
          u0 = self.u[i].T 
          
          u1 = u0 - np.linalg.inv(np.eye(len(self.dydt.u0)) - self.dt * self.dydt.df(self.time[i+1],u0)).dot(u0 - self.dt * self.dydt.f(self.time[i+1], u0).T - self.u[i].T )
          
          error = np.array([1.0])
          iters = 0
          while True:
              try:  
                 u0 = u1.T
                 u1 = u0 - np.linalg.inv(np.eye(len(self.dydt.u0)) - self.dt * self.dydt.df(self.time[i+1],u0)).dot(u0 - self.dt * self.dydt.f(self.time[i+1], u0).T - self.u[i].T )
                 
                 iters += 1
                 error = np.abs(u1-u0)
                 
                 if np.sum(np.abs(error)) <= toll:
                     break
                 if iters >= 1000:
                    raise ValueError('Too many iterations!') 
              except ValueError as ex:
                 print('Error occurred in Implicit-NR Euler Solvers: %s' %ex.args)
                 return None , None
          self.u[i+1] = u1.T
#-------------------------------------------      
        print('... Done!')
        Backward.solved = True 
        if Backward.solved and self.dydt.solution[0] != 0:
          self.dydt.evalError(self.u, 'Implicit NR-Euler', prints = False)


        if self.file != None:
            super().write2file()

        if self.save:
          return self.time,self.u

    def plot(self):
        '''
         Check if the problem is solved (if the solve method was already called) 
         and if it is plot the numerical solution[0] of the IV differential problem
        '''
        if Backward.solved:
            super().plot('ODE solution[0] using (Implicit) Backward Euler', 'time [s]', 'y(t)')    
        else:
            print("Unsolved problem, call `solve` method before")

#------------------------------------------------------------------------------------------------------
#   stepper

    @classmethod 
    def step(cls, func, t : 'np.float' , u : 'np.array', dt ):
        toll = 10e-12
        def f(t,u):
            return  np.array([function(t,u) for function in func])     
        def Df(t,u):
            eps = 10e-6
            return ((f(t,u)+eps) - f(t,u))/eps
        
        def df(t, u, **params):
            eps = 1e-12
            J = np.zeros([len(u), len(u)], dtype = np.float)
            for i in range(len(u)):
                u1 = u.copy()
                u2 = u.copy()
                u1[i] += eps
                u2[i] -= eps
                f1 = f(t, u1, **params)
                f2 = f(t, u2, **params)
                J[ : , i] = (f1 - f2) / (2 * eps)
            return J
    



        u0 = u.T
        u1 = u0 - np.linalg.inv(np.eye(len(u0)) - dt* df(t+dt,u0)).dot(u0 - dt* f(t+dt,u0).T - u.T )
        iters = 0
        while True:
           try:  
              u0 = u1.T
              u1 = u0 - np.linalg.inv(np.eye(len(u0)) - dt* df(t+dt,u0)).dot(u0 - dt* f(t+dt,u0).T - u.T )
              iters+=1
              if np.abs(u1-u0) <= toll:
                break
              if iters >= 1000:
                 raise ValueError('Too many Iterations!')
           except ValueError as ex :
               print('Error occurred in Implicit-NR Euler Stepper: %s', ex.args)
               return None
           #else
        return u1.T         

#------------------------------------------------------------------------------------------------------


