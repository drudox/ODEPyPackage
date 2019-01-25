# -*- coding : utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from ... Rhs.rhs import Rhs
from scipy.optimize import fsolve , newton_krylov ,newton
from .. RungeKutta import rungekutta #import RungeKutta
from .. RungeKutta.rungekutta import RungeKutta

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

class RadauIIA(RungeKutta):
    ''' 
        Solver for stiff and non-stiff differential equations
        solve the Problem using 3th order accuracy Implicit RadauIIA method
        
    '''
    solved = False
        
    def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            self.m = len(dydt.func)
            
            super().__init__(dydt)

    def implicit_equations(self, variables, ti, ui):
        k1, k2 = variables[:self.m], variables[self.m:]
        f1 = -k1 + ui + self.dt*(5./12 * self.f(ti + (1./3) * self.dt, k1) - 1/12. * self.f(ti + self.dt, k2 ) )
        f2 = -k2 + ui + self.dt*(3./4  * self.f(ti + (1./3) * self.dt, k1) + 1/4.  * self.f(ti + self.dt, k2 ) )  
        return np.hstack([f1,f2]).ravel()
  
    def solve(self): 
        s = 2 # s =2, the number of stadies of the num method
        self.time , self.u = self.dydt.createArray()
        print("Running Implicit RADAU-IIA order") 
        for i in range(len(self.time) - 1):
            res = fsolve(self.implicit_equations, np.random.rand(self.m * s), args=(self.time[i], self.u[i])) # s =2, the number of stadies of the num method
            k1, k2 = res[:self.m], res[self.m:]
            self.u[i + 1] = self.u[i] + self.dt * ( 3/4. * self.f(self.time[i] + (1./3)*self.dt ,k1) + 1/4. *self.f(self.time[i]+self.dt, k2))
        print('The system is solved.')
        print('... Done!')
        return self.time,self.u
    
        RadauIIA.solved = True
            

        if RadauIIA.solved and self.dydt.solution[0] != 0:
            self.dydt.evalError(self.u, 'Imp RK-4th', prints = False)
               
        if self.file != None:
            super().write2file()


        if self.save:
                return self.time,self.u

    def plot(self):
            if RadauIIA.solved:
               super().plot('ODE Solution using Radau 3th order','time [s]', 'y(t)')
            else:
               print("Unsolved problem, call `solve` method before")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

class RadauIIA5th(RungeKutta):
    ''' 
        Solver for stiff and non-stiff differential equations
        solve the Problem using 5th order accuracy Implicit RadauIIA 
        
    '''
    solved = False
        
    def __init__(self, dydt : Rhs , filename: str=None , save: bool=True ):
            self.file = filename 
            self.save = save
            self.m = len(dydt.func)
            
            super().__init__(dydt)

    def implicit_equations(self, variables, ti, ui):
        k1, k2 , k3 = variables[:self.m], variables[self.m:2*self.m] ,variables[2*self.m:]
        c1 = 2./5 - np.sqrt(6)/10.
        c2 = 2./5 + np.sqrt(6)/10.
        c3 = 1.
        a11 = 11./45 - 7*np.sqrt(6)/360.
        a12 = 37./225 - 169*np.sqrt(6)/1800.
        a13 = -2./225 + np.sqrt(6)/75
        a21 = 37./225 + 169*np.sqrt(6)/1800.
        a22 = 11./45  +   7*np.sqrt(6)/360.
        a23 = -2./225. - np.sqrt(6)/75.
        a31 = 4/9. - np.sqrt(6)/36. 
        a32 = 4/9. + np.sqrt(6)/36. 
        a33 = 1/9.

        f1 = -k1 + ui + self.dt*( a11 * self.f(ti+c1*self.dt , k1) + a12 * self.f(ti+c2*self.dt,k2) + a13*self.f(ti+c3*self.dt, k3))
        f2 = -k2 + ui + self.dt*( a21 * self.f(ti+c1*self.dt , k1) + a22 * self.f(ti+c2*self.dt,k2) + a23*self.f(ti+c3*self.dt, k3))  
        f3 = -k3 + ui + self.dt*( a31 * self.f(ti+c1*self.dt , k1) + a32 * self.f(ti+c2*self.dt,k2) + a33*self.f(ti+c3*self.dt, k3))
        return np.hstack([f1,f2,f3]).ravel()
  
    def solve(self): 
        s = 3 # s =3, the number of stadies of the num method
        self.time , self.u = self.dydt.createArray()
        print("Running Implicit RADAU-IIA 5th order") 
        for i in range(len(self.time) - 1):
            res = fsolve(self.implicit_equations, np.random.rand(self.m * s), args=(self.time[i], self.u[i]), xtol=10e-6) # s =2, the number of stadies of the num method
            k1,k2,k3 = res[:self.m], res[self.m:2*self.m], res[2*self.m:]
            b1 = 4./9 - np.sqrt(6)/36. 
            b2 = 4./9 + np.sqrt(6)/36.
            b3 = 1/9.
            c1 = 2./5 - np.sqrt(6)/10.
            c2 = 2./5 + np.sqrt(6)/10.
            c3 = 1.

            self.u[i + 1] = self.u[i] + self.dt * ( b1 * self.f(self.time[i] + c1*self.dt , k1) + b2 *self.f(self.time[i] + c2*self.dt, k2) + b3* self.f(self.time[i]+c3*self.dt, k3) )
        print('The system is solved.')
        print('... Done!')
        return self.time,self.u
    
        RadauIIA5th.solved = True
            

        if RadauIIA5th.solved and self.dydt.solution[0] != 0:
            self.dydt.evalError(self.u, 'Imp RadauIIA 5th', prints = False)
               
        if self.file != None:
            super().write2file()


        if self.save:
                return self.time,self.u

    def plot(self):
            if RadauIIA5th.solved:
               super().plot('ODE Solution using RadauIIA 5th order','time [s]', 'y(t)')
            else:
               print("Unsolved problem, call `solve` method before")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

 
