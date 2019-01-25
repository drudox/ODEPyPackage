# -*- coding : utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 
from . import euler
from ... Rhs.rhs import Rhs 
import scipy.optimize as opt
from scipy.optimize import fsolve , newton_krylov

class BWDEuler(euler.Euler):

    solved = False 

    def __init__(self, dydt : Rhs , filename : str= '', save : bool = True):
        self.save = save
        self.file = filename 
        
        super().__init__(dydt)
    
    def solve(self):

        self.time , self.u = self.dydt.createArray()

        print ('running Backward Euler ...')

        for i in range(len(self.time)-1):
            
            func = lambda up1 : up1 - self.u[i]-self.dt*self.f(self.time[i+1], up1) 
            self.u[i+1] = newton_krylov(func , self.u[i] , f_tol = 10e-8) 
                
        BWDEuler.solved = True  
        if BWDEuler.solved and self.dydt.solution[0] != 0:
           self.dydt.evalError(self.u, 'BWDEuler', prints=False)

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

        





