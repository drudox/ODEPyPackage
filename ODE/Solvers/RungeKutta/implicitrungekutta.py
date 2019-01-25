# -*- coding : utf-8 -*- 

import numpy as np
import matplotlib.pyplot as plt
from ... Rhs.rhs import Rhs
from scipy.optimize import fsolve , newton_krylov
from . rungekutta import RungeKutta


##-------------------------------------------------------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------------------------------------------------------
##                                         IMPLICIT RUNGE KUTTA METHODS  (STIFF SysODE SOLVER)
##-------------------------------------------------------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------------------------------------------------------



class ImpRK3(RungeKutta):

    solved = False
        
    def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            super().__init__(dydt)

    def solve(self):
        '''
              perform Implicit RK2 scheme 
        '''
        self.time , self.u = self.dydt.createArray()
        print("Running Implicit RK 2nd order") 
        
        for i in range(len(self.time)-1):

            func = lambda up1 : up1 - self.u[i] - 0.5*self.dt *(self.f(self.time[i],self.u[i]) + self.f(self.time[i+1] , up1) ) - 0.5* self.dt**2 *(self.dydt.Df(self.time[i],self.u[i]) - self.dydt.Df(self.time[i+1],up1) )
            
            self.u[i+1] = fsolve(func , self.u[i] , xtol = 10e-8) 
            #self.u[i+1] = newton_krylov(func , self.u[i] , f_tol = 10e-8) 

        print('... Done!')
        ImpRK3.solved = True
            

        if ImpRK3.solved and self.dydt.solution[0] != 0:
            self.dydt.evalError(self.u, 'Imp RK-2nd', prints = False)
               
        if self.file != None:
            super().write2file()


        if self.save:
                return self.time,self.u

    def plot(self):
            if ImpRK3.solved:
               super().plot('ODE Solution using Runge Kutta 2nd order','time [s]', 'y(t)')
            else:
               print("Unsolved problem, call `solve` method before")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

class ImpRK2(RungeKutta):

    solved = False
        
    def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            super().__init__(dydt)

    def solve(self):
        '''
              perform Implicit RK2 scheme 
        '''
        self.time , self.u = self.dydt.createArray()
        print("Running Implicit RK 2nd order") 
        
        for i in range(len(self.time)-1):

            func = lambda up1 : up1 - self.u[i] - self.dt/6 * ( 4.0* self.f(self.time[i], self.u[i]) + 2.* self.f(self.time[i+1],up1) + self.dt*self.dydt.Df(self.time[i],self.u[i]))
            self.u[i+1] = self.u[i+1] = fsolve(func , self.u[i] , xtol = 10e-8) 

        print('... Done!')
        ImpRK2.solved = True
            

        if ImpRK2.solved and self.dydt.solution[0] != 0:
            self.dydt.evalError(self.u, 'Imp RK-2nd', prints = False)
               
        if self.file != None:
            super().write2file()


        if self.save:
                return self.time,self.u

    def plot(self):
            if ImpRK2.solved:
               super().plot('ODE Solution using Runge Kutta 2nd order','time [s]', 'y(t)')
            else:
               print("Unsolved problem, call `solve` method before")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

class ImpRK4(RungeKutta):

    solved = False
        
    def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            self.m = len(dydt.func)
            
            super().__init__(dydt)

    def implicit_equations(self, variables, ti, ui):
        k1, k2 = variables[:self.m], variables[self.m:]
        f1 = -k1 + self.f(ti + (0.5 + np.sqrt(3) / 6) * self.dt, ui + 0.25 * self.dt * k1 + (0.25 + np.sqrt(3) / 6) * self.dt * k2)
        f2 = -k2 + self.f(ti + (0.5 - np.sqrt(3) / 6) * self.dt, ui + (0.25 - np.sqrt(3) / 6) * self.dt * k1 + 0.25 * self.dt * k2)
        return np.hstack([f1,f2]).ravel()
  
    def solve(self): 
        s = 2 # s =2, the number of stadies of the num method
        self.time , self.u = self.dydt.createArray()
        print("Running Implicit RK 4th order") 
        for i in range(len(self.time) - 1):
            res = fsolve(self.implicit_equations, np.random.rand(self.m * s), args=(self.time[i], self.u[i])) # s =2, the number of stadies of the num method
            k1, k2 = res[:self.m], res[self.m:]
            self.u[i+1] = self.u[i] + self.dt / 2 * (k1 + k2)
        print('The system is solved.')
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

       
