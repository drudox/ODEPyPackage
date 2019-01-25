'''
    module containing the adams method (Implicit and Explicit)
    @Marco Ghaini phD Glasgow Jul 2018
'''

from abc import abstractmethod
from .. ABSolver.sode import Sode
from .. RungeKutta import rungekutta
import numpy as np
import types
import matplotlib.pyplot as plt
from ... Rhs.rhs import Rhs
from .. MultiStep.multistep import MultiStep
from .. Euler import euler
#----------------------------------------------------------------------------------------------------------

class AdamsMethods(MultiStep):
    '''
        Base class for Adams-Bashfort / Adams-Moulton Solver

    '''
    @abstractmethod 
    def solve(self):
        pass

#----------------------------------------------------------------------------------------------------------
#   
#
#
#----------------------------------------------------------------------------------------------------------

class AdamsBashforth(AdamsMethods):
    '''     
       Implement several Adam-Bashforth multi-step explicit-scheme  
       Suitable for System of ODE (IVP) and for single differential 
       equations
    '''
    def __init__(self, dydt : Rhs, filename :str = None , save : bool = True):
        self._2nd = AdamsBashforth._2nd(dydt,filename,save)
        self._3rd = AdamsBashforth._3rd(dydt,filename,save)
        self._4th = AdamsBashforth._4th(dydt,filename,save)
        self._5th = AdamsBashforth._5th(dydt,filename,save)

    #@abstractmethod 
    def solve(self):
        pass
#--------------------------------------- 2nd order solver --------------------------------------------   
    class _2nd(AdamsMethods):
           startup = True
           solved  = False 
           um1 = 0

           def __init__(self, dydt : Rhs, filename :str = None , save : bool = True ):
                self.file = filename
                self.save = save
                super().__init__(dydt)
 
           def solve(self):
                print("Running Adams-Bashforth 2-step ....")
                self.time , self.u = self.dydt.createArray()
                               
                # RK2 start up the method 
                k1 = self.dydt.f(self.time[0],self.u[0])
                k2 = self.dydt.f(self.time[0]+ self.dt , self.u[0] + self.dt*k1 )
                self.u[1] = self.u[0] + self.dt/2. * (k1 + k2) 

                for i in range(1,len(self.time)-1):
                    self.u[i+1] = self.u[i] + self.dt/2.*(3.*self.dydt.f(self.time[i],self.u[i])-self.dydt.f(self.time[i-1],self.u[i-1]))
                
                print('... Done!')
                AdamsBashforth._2nd.solved = True    
                if AdamsBashforth._2nd.solved and self.dydt.solution != None:
                    self.dydt.evalError(self.u, 'AdamsBashforth 2nd', prints = False)
    
               
                if self.file != None:
                    super().write2file()
            
                if self.save:
                    return self.time,self.u

                def plot(self):
                    if AdamsBashforth._2nd.solved:
                        super().plot('Sys ODE solution using Adams-Bashforth 2nd order','time [s]','y(t)')
                    else:
                        print("Unsolved problem, call `solve` method before")
 
#-----------------------------------------------------------------------------------------------------------

           @classmethod
           def step(cls, func , t : np.float , u : np.float, dt ):
        
                def f(ti,ui):
                    return  np.array([function(ti,ui) for function in func])     
                            
                if cls.startup:
                    cls.um1 = u.copy()    # u[i-1]
                    unext = rungekutta.RK2.Heun.step(func,t,cls.um1,dt)  #cls.um1 + dt/2.*(k1+k2)
                    cls.startup = False
                else:
                    unext = u + dt/ 2.* (3.* f(t,u) - f(t-dt, cls.um1))
                    cls.um1 = u
                return unext 

###---------------------------------------------------------------------------------------------------------
###                               Adams Bashforth   3 step          
###---------------------------------------------------------------------------------------------------------
    class _3rd(AdamsMethods):
           first_startup = True
           second_startup = True
           solved  = False 
           um1 = 0.
           um2 = 0.

           def __init__(self, dydt : Rhs, filename :str = None , save : bool = True ):
                self.file = filename
                self.save = save
                super().__init__(dydt)
 
           def solve(self):
                print("Running Adams-Bashforth 3-step ....")
                self.time , self.u = self.dydt.createArray()
                
                # RK2 start up the method 
                k1 = self.dydt.f(self.time[0],self.u[0])
                k2 = self.dydt.f(self.time[0] + self.dt , self.u[0] +self.dt * k1)
                self.u[1] = self.u[0] + self.dt/2. * (k1 + k2)
                
                # second point 
                k1   = self.dydt.f(self.time[1],self.u[1])
                k2   = self.dydt.f(self.time[1] +self.dt , self.u[1]+ self.dt * k1)
                self.u[2] = self.u[1] + self.dt/2. * (k1 + k2)
                
                for i in range(2,len(self.time)-1):
                    self.u[i+1] = self.u[i] + self.dt/12.*(23.*self.dydt.f(self.time[i],self.u[i]) \
                                                          -16.*self.dydt.f(self.time[i-1],self.u[i-1]) \
                                                          + 5.*self.dydt.f(self.time[i-2],self.u[i-2]) )
                
                print('...Done!')
                AdamsBashforth._3rd.solved = True
                if AdamsBashforth._3rd.solved and self.dydt.solution != None:
                    self.dydt.evalError(self.u, 'AdamsBashforth 3rd', prints = False)
    
 
                if self.file != None:
                    super().write2file()
            
                if self.save:
                    return self.time,self.u
                
                def plot(self):
                    if AdamsBashforth._3rd.solved:
                        super().plot('Sys ODE solution using Adams-Bashforth 3rd order','time [s]','y(t)')
                    else:
                        print("Unsolved problem, call `solve` method before")
#-----------------------------------------------------------------------------------------------------------
#                                 Adams-Bashfort 3th order  Stepper 
           @classmethod
           def step(cls, func , t :np.float, u :np.array, dt):
                
                def f(ti,ui):
                    return  np.array([function(ti,ui) for function in func])  
                
                if cls.first_startup:
                    cls.um2 = u.copy()
                    unext = rungekutta.RK2.Heun.step(func,t,cls.um2,dt)  #cls.um1 + dt/2.*(k1+k2)
                    cls.first_startup = False
                    t+=dt
                elif cls.second_startup:
                    cls.um1 = u.copy()
                    unext = rungekutta.RK2.Heun.step(func,t,cls.um1,dt)  #cls.um1 + dt/2.*(k1+k2)
                    cls.second_startup = False
                    t+=dt
                else:  # compute Adams Bashford    
                    unext = u + dt/12.*(23*f(t,u) -16*f(t-dt,cls.um1) + 5*f(t-dt-dt,cls.um2))
                    cls.um2 = cls.um1.copy()
                    cls.um1 = u.copy()
                return unext
###---------------------------------------------------------------------------------------------------------
###                                     Adams Bashforth   4 step          
###---------------------------------------------------------------------------------------------------------
    
    class _4th(AdamsMethods):
           '''
                  4th order 4-step explicit Adams-Bashforth 
           '''
           um1 , um2 , um3 = np.array([0.]), np.array([0.]), np.array([0.])
           first_startup, second_startup, third_startup = True , True, True
           solved = False
           def __init__(self, dydt : Rhs , filename :str = ' ', save :bool = True):
               self.file = filename
               self.save = save
               super().__init__(dydt)
           
           def solve(self):
               
               print('Running Adams-Bashforth 4-step ....')
               
               self.time, self.u = self.dydt.createArray()
               
               # start-up using RK4 high order  
               for i in range(3):    
                   k1 = self.dydt.f(self.time[i], self.u[i])
                   k2 = self.dydt.f(self.time[i] + self.dt/2., self.u[i] + self.dt/2. *k1 )
                   k3 = self.dydt.f(self.time[i] + self.dt/2., self.u[i] + self.dt/2. *k2 )
                   k4 = self.dydt.f(self.time[i] + self.dt   , self.u[i] + self.dt    *k3 )       
                   
                   self.u[i+1] = self.u[i] + self.dt/6.* (k1 + 2*k2 + 2*k3 + k4)      

               for i in range(3,len(self.time)-1):    
                   #
                   self.u[i+1] = self.u[i] + self.dt/24. * (55. * self.dydt.f(self.time[i], self.u[i]) \
                                                          - 59. * self.dydt.f(self.time[i-1], self.u[i-1])  \
                                                          + 37. * self.dydt.f(self.time[i-2], self.u[i-2]) \
                                                          -  9  * self.dydt.f(self.time[i-3], self.u[i-3]) )
                   
               print('... Done!')
               AdamsBashforth._4th.solved = True
               if AdamsBashforth._4th.solved and self.dydt.solution != None:
                    self.dydt.evalError(self.u, 'AdamsBashforth 4th', prints = False)
    
               if self.file != None:
                  super().write2file() 
                  
               if self.save:
                  return self.time,self.u
                   
               def plot(self):
                  if AdamsBashforth._4th.solved:
                     super().plot('Sys ODE solution using Adams-Bashforth 4th order','time [s]','y(t)')
                  else:
                     print("Unsolved problem, call `solve` method before")
#----------------------------------------------------------------------------------------------------------
           @classmethod
           def step(cls, func , t : np.float , u : np.float , dt):
             
                def f(ti,ui):
                    return  np.array([function(ti,ui) for function in func])  
            
                if cls.first_startup:
                    cls.um3 = u.copy()
                    unext = rungekutta.RK4.step(func,t,cls.um3,dt) 
                    t += dt
                    cls.first_startup = False
                elif cls.second_startup:
                    cls.um2 = u.copy()
                    unext = rungekutta.RK4.step(func,t,cls.um2,dt) 
                    t+= dt
                    cls.second_startup = False
                elif cls.third_startup:  
                    cls.um1 = u.copy()
                    unext = rungekutta.RK4.step(func,t,cls.um1,dt)
                    t += dt 
                    cls.third_startup = False
                else:  # compute AB 4th order 
                    unext = u + dt/24.* ( 55.*f(t,u) - 59.*f(t-dt,cls.um1) + 37.*f(t-(2*dt),cls.um2) \
                                                                           - 9.*f(t-(3*dt),cls.um3 )) 
                    cls.um3 = cls.um2.copy()
                    cls.um2 = cls.um1.copy()
                    cls.um1 = u.copy() 
                return unext   
#-----------------------------------------------------------------------------------------------------------
#  up = u + dt/24. * (55.*f(t,u) - 59*f(t-dt,cls.um1) + 37.*f(t-(2*dt), cls.um2) + f(t-(3*dt),cls.um3) ) 
                   
#
#Adams Bashforth 5-step (5th order)
#-----------------------------------------------------------------------------------------------------------
    class _5th(AdamsMethods):
        
        um4,um3,um2,um1 = ([np.array]),([np.array]),([np.array]),([np.array])
        fst_startup, snd_startup, trd_startup, fth_startup = True,True,True,True
        solved = False

        def __init__(self, dydt : Rhs , filename :str=' ', save: bool=True ):
            self.file = filename
            self.save = save
            super().__init__(dydt)

        def solve(self):
            '''

            '''

            self.time,self.u = self.dydt.createArray()
            print('Running Adams-Bashforth 5-step ....')

            # compute the first fourth points
            for i in range(4):
                k1 = self.dydt.f(self.time[i] , self.u[i])
                k2 = self.dydt.f(self.time[i] + self.dt/2. , self.u[i] +self.dt/2.*k1)
                k3 = self.dydt.f(self.time[i] + self.dt/2. , self.u[i] +self.dt/2.*k2)
                k4 = self.dydt.f(self.time[i] + self.dt    , self.u[i] +self.dt   *k3)

                self.u[i+1] = self.u[i]+ self.dt/6. * (k1 + 2* k2 + 2* k3 + k4) 
            
            for i in range(4,len(self.time)-1):
                self.u[i+1] = self.u[i] + self.dt * (1901./720. * self.dydt.f(self.time[i], self.u[i]) \
                                                    -1387./360. * self.dydt.f(self.time[i-1], self.u[i-1]) \
                                                    + 109./30.  * self.dydt.f(self.time[i-2], self.u[i-2])\
                                                    - 637./360. * self.dydt.f(self.time[i-3], self.u[i-3])\
                                                    + 251./720. * self.dydt.f(self.time[i-4], self.u[i-4]) )
            print('... Done!')
            AdamsBashforth._5th.solved = True      
            if AdamsBashforth._5th.solved and self.dydt.solution != None:
                 self.dydt.evalError(self.u, 'AdamsBashforth 5th', prints = False)
    

            print(self.file != None)                            
            if self.file != None:
                super().write2file() 
                  
            if self.save:
                return self.time,self.u
                   
            def plot(self):
                if AdamsBashforth._5th.solved:
                    super().plot('Sys ODE solution using Adams-Bashforth 5th order','time [s]','y(t)')
                else:
                    print("Unsolved problem, call `solve` method before")

#----------------------------------------------------------------------------------------------------------- 
#   stepper
    
        @classmethod
        def step(cls, func , t, u, dt):
            
            def f(ti,ui):
                return np.array([function(ti,ui) for function in func])
            
            if cls.fst_startup:
                    cls.um4 = u.copy()
                    unext   = rungekutta.RK4.step(func,t,cls.um4,dt)
                    t += dt
                    cls.fst_startup = False
            elif cls.snd_startup:
                    cls.um3 = u.copy()
                    unext   = rungekutta.RK4.step(func,t,cls.um3,dt)
                    t+=dt
                    cls.snd_startup = False
            elif cls.trd_startup:
                    cls.um2 = u.copy()
                    unext   = rungekutta.RK4.step(func,t,cls.um2,dt)
                    t += dt
                    cls.trd_startup = False
            elif cls.fth_startup:
                    cls.um1 = u.copy()
                    unext   = rungekutta.RK4.step(func,t,cls.um1,dt)
                    t += dt
                    cls.fth_startup = False
            else :
                    unext = u + dt * (1901./720. * f(t,u) - 1387./360. * f(t-dt, cls.um1) \
                                      + 109./30. * f(t-(2*dt),cls.um2) - 637./360. * f(t-(3*dt), cls.um3)\
                                      + 251./720. * f(t-(4*dt),cls.um4) )
                    
                    cls.um4 = cls.um3
                    cls.um3 = cls.um2
                    cls.um2 = cls.um1
                    cls.um1 = u
                
            return unext    



############################################################################################################
###------------------------------------------------------------------------------------------------------###
###                                                                                                      ###       
###                               ADAMS BASHFORTH-MOULTON SEMI-IMPLICIT SCHEME                           ###   
###                                                                                                      ###
###------------------------------------------------------------------------------------------------------###
############################################################################################################

class AdamsMoulton(AdamsMethods):
    '''

    '''
    def __init__(self, dydt, filename :str = None , save : bool = True):
        self._2nd = AdamsMoulton._2nd(dydt,filename,save)
        self._3rd = AdamsMoulton._3rd(dydt,filename,save)
        self._4th = AdamsMoulton._4th(dydt,filename,save)
        self._5th = AdamsMoulton._5th(dydt,filename,save)
    
    def solve(self):
        pass
#--------------------------------------- 2nd order solver --------------------------------------------   
    class _2nd(AdamsMethods):
        ''' 
            Perform the predictor - corrector solutions (Semi implicit Predictor-Corrector)
            of a System of Ordinary Differential Equations
        '''
        um1, toll, startup, solved = 0., 1E-12, True, False

        def __init__(self, dydt : Rhs, filename : str=None , save : bool = True):
            
            self.file = filename
            self.save = save
            super().__init__(dydt)
        
        def solve(self):
                print("Running Adams-Moulton 2nd order (Implicit) ....")
                self.toll = 1E-12
                self.time , self.u = self.dydt.createArray()
                
                up = self.u.copy()
                # Start-Up The Predictor (AB2 method)
                k1 = self.dydt.f(self.time[0],self.u[0])
                k2 = self.dydt.f(self.time[0]+ self.dt , self.u[0] + self.dt*k1 )
                self.u[1] = self.u[0] + self.dt/2. * (k1 + k2) 

                for i in range(1,len(self.time)-1):
                    
                    #--------------------    PREDICTOR STEP    -----------------
                    up[i] = self.u[i] + self.dt/2.*(3.*self.dydt.f(self.time[i],self.u[i])-self.dydt.f(self.time[i-1],self.u[i-1]))
                    
                    #------------------ CORRECTOR (IMPLICIT A-M) ---------------- 
                    fpred = self.dydt.f(self.time[i+1],up[i])   
                    error = 1.
                    fold  = fpred
                    
                    iter = 0
                    while True:
                        uold = self.u[i] + self.dt/2.*(fold + self.dydt.f(self.time[i],self.u[i]))   

                        fold = self.dydt.f(self.time[i+1],uold)

                        unew = self.u[i] + self.dt/2.*(fold + self.dydt.f(self.time[i],self.u[i]))

                        fnew = self.dydt.f(self.time[i+1],unew)

                        error = np.abs(unew-uold)
                        fold  = fnew 
                        iter += 1
                        if error[0] <= self.toll: 
                            break
                    
                    self.u[i+1] = unew
                
                print('... Done!')
                
                AdamsMoulton._2nd.solved = True    
                
                if AdamsMoulton._2nd.solved and self.dydt.solution != None:
                    self.dydt.evalError(self.u, 'AdamsMoulton 2nd', prints = False)
    
                if self.file != None:
                    super().write2file()
            
                if self.save:
                    return self.time,self.u

        def plot(self):
            if AdamsMoulton._2nd.solved:
                super().plot('Sys ODE solution using Adams-Bashforth 2nd order','time [s]','y(t)')
            else:
                print("Unsolved problem, call `solve` method before")
#----------------------------------------------------------------------------------------------------------- 

        @classmethod
        def step(cls, func , t : np.float , u : np.float, dt ):
            
            def f(ti,ui):
                return  np.array([function(ti,ui) for function in func])     
               
            if cls.startup:
                cls.um1 = u.copy()
                unext = rungekutta.RK2.Heun.step(func,t,cls.um1,dt) #u + dt/2*(k1+k2)
                cls.startup = False
            else:
                # PREDICTOR
                up = u + dt/2.*(3* f(t,u) - f(t-dt,cls.um1))
                cls.um1 = u.copy()
                # CORRECTOR 
                fpred = f(t+dt,up)
                error = 1.
                fold  = fpred
                iter = 0
                while True:
                    uold = u+ dt/2.*(fold+ f(t,u))
                    fold = f(t+dt,uold)
                    unew = u+ dt/2.*(fold+ f(t,u))
                    fnew = f(t+dt, unew)
                    error= np.abs(unew-uold)
                    fold = fnew
                    iter +=1 
                    if np.sum(np.abs(error)) <= cls.toll:
                       break
                unext = unew
            return unext    

#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#                          Adams Moulton 2 Step (3th order implicit scheme) 
#                  @Compute the Predictor - Implicit Corrector method to solve a System of ODE
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
    
    class _3rd(AdamsMethods):
        
        first_startup, second_startup, solved = True, True , False 
        um1,um2 = np.array([0]) , np.array([0.])
        toll = 1E-12
        def __init__(self, dydt : Rhs, filename = None , save : bool = True ):
             self.file = filename
             self.save = save
             super().__init__(dydt)
 
        def solve(self):

             print("Running Adams-Moulton (Implicit) 3th-ord ....")
             self.time , self.u = self.dydt.createArray()
                
                # RK2 start up the method 
             k1 = self.dydt.f(self.time[0],self.u[0])
             k2 = self.dydt.f(self.time[0] + self.dt , self.u[0] +self.dt * k1)
             self.u[1] = self.u[0] + self.dt/2. * (k1 + k2)
                
                # second point 
             k1   = self.dydt.f(self.time[1],self.u[1])
             k2   = self.dydt.f(self.time[1] +self.dt , self.u[1]+ self.dt * k1)
             self.u[2] = self.u[1] + self.dt/2. * (k1 + k2)
              
              # MAiN LOOP 
             for i in range(2,len(self.time)-1):
                  '''
                        main loop, here is compute the Adams Bashforth 3th ord. predictor
                        with the Implicit scheme Adams Moulton 3th ord. Corrector 
                        iterating over the corrector step
                  '''  
                  #-----------------------------    PREDICTOR STEP    ----------------------------------
                  uold  = self.u[i] + self.dt/12.*(23.*self.dydt.f(self.time[i],self.u[i]) \
                                                  -16.*self.dydt.f(self.time[i-1],self.u[i-1]) \
                                                  + 5.*self.dydt.f(self.time[i-2],self.u[i-2]) )
                  
                  fpred = self.dydt.f(self.time[i+1],uold)
                  error = 1.0  
                  fold  = fpred
                  
                  #----------------------------  CORRECTOR IMPLICIT STEP    -----------------------------
                  while True:
                    
                      uold = self.u[i] + self.dt/12. * ( 5* fold \
                                                       + 8* self.dydt.f(self.time[i], self.u[i]) \
                                                       - 1* self.dydt.f(self.time[i-1],self.u[i-1]) )  
                    
                      fold = self.dydt.f(self.time[i+1],uold)
                    
                      unew = self.u[i] + self.dt/12. * ( 5* fold \
                                                       + 8* self.dydt.f(self.time[i], self.u[i]) \
                                                       - 1* self.dydt.f(self.time[i-1],self.u[i-1]) )  
                      
                      fnew = self.dydt.f(self.time[i+1], unew)
                      error = abs(unew - uold)
                      fold = fnew  
                      if np.sum(np.abs(error)) <= self.toll:
                          break

                  self.u[i+1] = unew       
             
             #-----
              
             print('...Done!')
             AdamsMoulton._3rd.solved = True
             if AdamsMoulton._3rd.solved and self.dydt.solution != None:
                    self.dydt.evalError(self.u, 'AdamsMoulton 3rd', prints = False, write=True)
    
 
             if self.file != None:
                 super().write2file()
            
             if self.save:
                 return self.time,self.u
                
             def plot(self):
                if AdamsMoulton._3rd.solved:
                    super().plot('Sys ODE solution using Adams-Bashforth 2nd order','time [s]','y(t)')
                else:
                    print("Unsolved problem, call `solve` method before")
#----------------------------------------------------------------------------------------------------------
#              STEPPER  suitable for ODE         
#----------------------------------------------------------------------------------------------------------        
        @classmethod
        def step(cls , func , t :np.array , u :np.array, dt):
            
                def f(ti,ui):
                    return  np.array([function(ti,ui) for function in func])  
                
                if cls.first_startup:
                    cls.um2 = u.copy()
                    unext = rungekutta.RK2.Heun.step(func,t,cls.um2,dt)  #cls.um1 + dt/2.*(k1+k2)
                    cls.first_startup = False
                    t+=dt
                elif cls.second_startup:
                    cls.um1 = u.copy()
                    unext = rungekutta.RK2.Heun.step(func,t,cls.um1,dt)  #cls.um1 + dt/2.*(k1+k2)
                    cls.second_startup = False
                    t+=dt
                else:  # compute Adams Moulton
                    up = u + dt/12.*(23* f(t,u) - 16* f(t-dt,cls.um1) + 5* f(t-dt-dt,cls.um2))
                    cls.um2 = cls.um1.copy()
                    # CORRECTOR 
                    fpred = f(t+dt,up)
                    error = 1.
                    fold  = fpred
                    iter = 0
                    while True:
                        uold = u + dt/12.*(5*fold + 8* f(t,u) -1*f(t-dt, cls.um1))
                        fold = f(t+dt,uold)
                        unew = u + dt/12.*(5*fold + 8* f(t,u) -1*f(t-dt, cls.um1))
                        fnew = f(t+dt, unew)
                        error= np.abs(unew-uold)
                        fold = fnew
                        iter +=1 
                        if np.sum(np.abs(error)) <= cls.toll:
                            break
                    cls.um2 = cls.um1
                    cls.um1 = u.copy()
                    unext = unew
                return unext    
##----------------------------------------------------------------------------------------------------------
#
#                                      ADAMS-MOULTON  4th order 
#   @author Marco Ghiani
#-----------------------------------------------------------------------------------------------------------
 
    class _4th(AdamsMethods):
        
        um3,um2,um1 = np.array([0.]), np.array([0.]), np.array([0.])
        first_startup, second_startup, third_startup = True , True, True  
        solved = False
        toll = 1E-12
        def __init__(self, dydt : Rhs, filename = None , save :bool= True):
                self.file  = filename
                self.save  = save
                super().__init__(dydt)

        def solve(self):
            
            self.time, self.u = self.dydt.createArray()
            print('Running Adams-Moulton (Implicit) 4th-ord ...')
            
            # compute the first 3 point using RK-4 
            for i in range(3):  
                k1 = self.dydt.f(self.time[i],self.u[i])
                k2 = self.dydt.f(self.time[i]+ self.dt/2. , self.u[i] + self.dt/2.*k1 )
                k3 = self.dydt.f(self.time[i]+ self.dt/2. , self.u[i] + self.dt/2.*k2 )
                k4 = self.dydt.f(self.time[i]+ self.dt    , self.u[i] + self.dt  * k3 )

                self.u[i+1] = self.u[i] + self.dt/6. * (k1 + 2.*k2 + 2.*k3 + k4)
            
            ##------- AB4 PREDICTOR - AM4 CORRECTOR ----------
            
            for i in range(3,len(self.time)-1):
                '''

                '''
                #------------------------------ PREDICTOR STEP ------------------------------
                up = self.u[i] + self.dt/24. * ( 55. * self.dydt.f(self.time[i], self.u[i]) \
                                               - 59. * self.dydt.f(self.time[i-1],self.u[i-1])\
                                               + 37. * self.dydt.f(self.time[i-2],self.u[i-2])\
                                               -  9. * self.dydt.f(self.time[i-3],self.u[i-3]))
                
                fpred = self.dydt.f(self.time[i+1], up) 
                error = 1.0
                fold  = fpred
                #-------------------------- CORRECTOR IMPLICIT STEP---------------------------
                iter = 0
                while True:
                    
                    uold = self.u[i] + self.dt/24. * ( 9. * fold   \
                                                    + 19. * self.dydt.f(self.time[i], self.u[i]) \
                                                    -  5. * self.dydt.f(self.time[i-1],self.u[i-1])\
                                                    +  1. * self.dydt.f(self.time[i-2],self.u[i-2]) )
                    
                    fold = self.dydt.f(self.time[i+1], uold)  
                    
                    unew = self.u[i] + self.dt/24. * ( 9. * fold   \
                                                    + 19. * self.dydt.f(self.time[i], self.u[i]) \
                                                    -  5. * self.dydt.f(self.time[i-1],self.u[i-1])\
                                                    +  1. * self.dydt.f(self.time[i-2],self.u[i-2]) )
                                        
                    fnew = self.dydt.f(self.time[i+1], unew)
                    fold = fnew
                    error = np.abs(unew-uold)
                    iter += 1
                    if np.sum(np.abs(error))<= self.toll:
                        break
                
                self.u[i+1] = unew
                
            #------
            print('...Done!')    
            AdamsMoulton._4th.solved = True
            if AdamsMoulton._4th.solved and self.dydt.solution != None:
                  self.dydt.evalError(self.u, 'AdamsMoulton 4th', prints = False)
    

            if self.file != None:
                 super().write2file()
            
            if self.save:
                return self.time,self.u
                
            def plot(self):
               if AdamsMoulton._4th.solved:
                  super().plot('Sys ODE solution using Adams-Moulton 4th order','time [s]','y(t)')
               else:
                  print("Unsolved problem, call `solve` method before")

#----------------------------------------------------------------------------------------------------------
#     stepper 4th order 
        @classmethod
        def step( cls , func , t : np.float , u : np.array , dt):
            
            def f(ti,ui):
                return np.array([function(ti,ui) for function in func])

            # compute the first 3 point (start-up predictor solver) 

            if cls.first_startup:
                cls.um3 = u.copy()
                unext   = rungekutta.RK4.step(func,t,cls.um3,dt)
                t+=dt
                cls.first_startup = False
            elif cls.second_startup:
                cls.um2 = u.copy() 
                unext   = rungekutta.RK4.step(func,t,cls.um2,dt)
                t +=dt
                cls.second_startup = False
            elif cls.third_startup:
                cls.um1 = u.copy()
                unext   = rungekutta.RK4.step(func,t,cls.um1,dt)
                t +=dt
                cls.third_startup = False
            else:
                up = u + dt/24.* (55. * f(t,u) - 59. * f(t-dt,cls.um1) + 37.*f(t-(2*dt),cls.um2) \
                                 - 9. * f(t-(3*dt),cls.um3))
                fpred = f(t+dt,up)  #
                error = 1.
                fold  = fpred
                iter = 0
                while True:
                    uold = u + dt/24. * (9.*fold + 19.*f(t,u)-5*f(t-dt,cls.um1)+1*f(t-(2*dt),cls.um2))
                    fold = f(t+dt, uold)
                    unew = u + dt/24. * (9.*fold + 19.*f(t,u)-5*f(t-dt,cls.um1)+1*f(t-(2*dt),cls.um2))
                    fnew = f(t+dt, unew)
                    error= np.abs(unew-uold)
                    if error <= cls.toll:
                        break
                cls.um3 = cls.um2
                cls.um2 = cls.um1
                cls.um1 = u
                unext = unew
            return unext 

#----------------------------------------------------------------------------------------------------------

    class _5th(AdamsMethods):
        '''
        '''
        um1,um2,um3,um4 = np.array([0]), np.array([0]), np.array([0]), np.array([0])
        fst_stup , snd_stup , trd_stup , fth_stup = True,True,True,True 
        toll = 1E-12

        def __init__(self, dydt : Rhs, filename = None, save : bool='True'):
            self.file = filename
            self.save = save
            super().__init__(dydt)
        
        def solve(self):
            
            self.time,self.u = self.dydt.createArray()
            print('Running Adams-Moulton (Implicit) 5th-ord')
            
            ## compute first fourth points using 4th order (classic) Runge Kutta

            for i in range(4):  
                k1 = self.dydt.f(self.time[i],self.u[i])
                k2 = self.dydt.f(self.time[i]+self.dt/2. , self.u[i]+ self.dt/2. *k1 )
                k3 = self.dydt.f(self.time[i]+self.dt/2. , self.u[i]+ self.dt/2. *k2 )
                k4 = self.dydt.f(self.time[i]+self.dt    , self.u[i]+ self.dt    *k3 )

                self.u[i+1] = self.u[i] + self.dt/6. * (k1 + 2*k2 + 2*k3 + k4)

            for i in range(4, len(self.time)-1):
                
                # PREDICTOR STEP AB5
                up = self.u[i] + self.dt * (1901./720. * self.dydt.f(self.time[i], self.u[i]) \
                                           -1387./360. * self.dydt.f(self.time[i-1], self.u[i-1]) \
                                           + 109./30.  * self.dydt.f(self.time[i-2], self.u[i-2])\
                                           - 637./360. * self.dydt.f(self.time[i-3], self.u[i-3])\
                                           + 251./720. * self.dydt.f(self.time[i-4], self.u[i-4]) )
                
                fpred = self.dydt.f(self.time[i+1], up)
                error = 1.
                fold  = fpred
                ## ------------------------- IMPLICIT CORRECTOR STEP --------------------------------
                iter = 0
                while True:
                    
                    uold = self.u[i] + self.dt/720. * ( 251. * fold  \
                                                      + 646. * self.dydt.f(self.time[i],self.u[i])\
                                                      - 264. * self.dydt.f(self.time[i-1],self.u[i-1])\
                                                      + 106. * self.dydt.f(self.time[i-2],self.u[i-2])\
                                                      -  19. * self.dydt.f(self.time[i-3],self.u[i-3]))
                    
                    fold = self.dydt.f(self.time[i+1],uold)
            
                    unew = self.u[i] + self.dt/720. * ( 251. * fold  \
                                                      + 646. * self.dydt.f(self.time[i],self.u[i])\
                                                      - 264. * self.dydt.f(self.time[i-1],self.u[i-1])\
                                                      + 106. * self.dydt.f(self.time[i-2],self.u[i-2])\
                                                      -  19. * self.dydt.f(self.time[i-3],self.u[i-3]))
                    
                    fnew = self.dydt.f(self.time[i+1],unew)
                    error = np.abs(unew-uold) 
                    fold = fnew
                    if np.sum(np.abs(error)) <= self.toll :
                        break
                self.u[i+1] = unew 
                
            print('... Done!')
            AdamsMoulton._5th.solved = True
            if AdamsMoulton._5th.solved and self.dydt.solution != None:
                  self.dydt.evalError(self.u, 'AdamsMoulton 5th', prints = False)
  
            if self.file != None:
                super().write2file()
           
            if self.save:
                return self.time,self.u
            
            def plot(self):
               if AdamsMoulton._5th.solved:
                  super().plot('Sys ODE solution using Adams-Moulton 5th order','time [s]','y(t)')
               else:
                  print("Unsolved problem, call `solve` method before")
 #----------------------------------------------------------------------------------------------------------


        @classmethod
        def step(cls, func, t : np.float, u : np.array , dt : np.float):
            
            def f(ti,ui):
                return np.array([function(ti,ui) for function in func ])

            # compute first 4 point using Rk4 stepper 
            if cls.fst_stup:
                cls.um4 = u.copy()
                unext   = rungekutta.RK4.step(func,t,cls.um4,dt)
                t+=dt
                cls.fst_stup = False
            elif cls.snd_stup:
                cls.um3 = u.copy()
                unext   = rungekutta.RK4.step(func,t,cls.um3,dt)
                t+=dt
                cls.snd_stup = False
            elif cls.trd_stup:
                cls.um2 = u.copy()
                unext   = rungekutta.RK4.step(func,t,cls.um2,dt)
                t += dt 
                cls.trd_stup = False
            elif cls.fth_stup:
                cls.um1 = u.copy()
                unext   = rungekutta.RK4.step(func,t,cls.um1,dt)
                t += dt
                cls.fth_stup = False
            else: # PREDICTOR ADAMS-BASHFORTH 5th order 
                up    = u + dt * (1901./720. * f(t,u) - 1387./360. * f(t-dt, cls.um1) \
                                  + 109./30. * f(t-(2*dt),cls.um2) - 637./360. * f(t-(3*dt), cls.um3)\
                                 + 251./720. * f(t-(4*dt),cls.um4) )
                fpred = f(t+dt,up)     
                error = 1.    
                fold  = fpred
                iter  = 0

                while True:
                    
                    uold = u + dt/720. * (251.* fold + 646.* f(t,u) - 264.*f(t-dt,cls.um1)  \
                                        + 106.* f(t-(2*dt),cls.um2) - 19. *f(t-(3*dt),cls.um3))
                    
                    fold = f(t+dt, uold)
                    
                    unew = u + dt/720. * (251.* fold + 646.* f(t,u) - 264.*f(t-dt,cls.um1)  \
                                        + 106.* f(t-(2*dt),cls.um2) - 19. *f(t-(3*dt),cls.um3))
                    fnew = f(t+dt, unew)
                    
                    error = np.abs(unew-uold)
                    fold  = fnew
                    if error <= cls.toll:
                        break


                cls.um4 = cls.um3
                cls.um3 = cls.um2
                cls.um2 = cls.um1
                cls.um1 = u
                unext   = unew
        
            return unext













