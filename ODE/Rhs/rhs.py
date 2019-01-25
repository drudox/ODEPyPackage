import types
import numpy as np
import matplotlib.pyplot as plt 
#import symengine as se                  # tools for Jacobian
#import numdifftools as nd

class Rhs:
   '''
      class to manage differential problem 
      contains the dy/dt function and the initial value 
      in order to close the ode problem
   '''
   solution = None
   analitical_called = False 
   def __init__(self, fnum : np.ndarray , t0: np.float, tf: np.float, y0 : np.array, n: int, fanal = np.array([0]), dt : np.float = None ):
      '''   
         Input : 
            - fnum : Function f(t,y(t)) = dy/dt (numeric)
            - 
            - Initial Time t0 
            - Final Time tf
            - Initial value y0(t0) 
      '''
      self.func = fnum
      self.solution = fanal  
      self.t0   = t0
      self.tf   = tf
      self.u0   = y0
      if dt != None:
        self.n = round((tf-t0)/dt)
      else:
        self.n = n  

   def createArray(self):
      '''   
         Create the Array time and f(time) 
         - the time array can be create now 
         - the solution array can be just initialize with the IV
      ''' 
      self.t = np.linspace(self.t0, self.tf, self.n+1 )
      self.u = np.array([self.u0  for i in range(self.n+1) ])
      
      
      return self.t,self.u   
  
   
   def exactsolution(self,xi,yi):
       return np.array([function(xi,yi)  for function in self.solution])   
    

   def f(self,ti,ui):
       return  np.array([function(ti,ui) for function in self.func])     
       

   def Df(self,ti,ui):
      eps = 10e-12
      return (self.f(ti,ui+eps) - self.f(ti,ui-eps) )/(2*eps)



   def df(self, t, u, **params):
      """
          Compute the Jacobian numerically
      """
      eps = 1e-12
      J = np.zeros([len(u), len(u)], dtype = np.float)

      for i in range(len(u)):
          u1 = u.copy()
          u2 = u.copy()

          u1[i] += eps
          u2[i] -= eps

          f1 = self.f(t, u1, **params)
          f2 = self.f(t, u2, **params)

          J[ : , i] = (f1 - f2) / (2 * eps)

      return J
    



#-----------------------------------------------------------------------------------------------------------
   #@staticmethod
   def analiticalSolution(self,save: 'if True (default) returm solution vectors'=False, show: 'plot soluition'=False ):
        
      if self.exactsolution != None: 
         self.x = np.linspace(self.t0,self.tf,self.n+1)
         self.y = np.array([self.u0 for i in range(self.n+1)])
         
         self.y = self.exactsolution(self.x,self.y).T 
         Rhs.analitical_called = True
         self.save = save
         self.show = show
         
         if self.save:
            return self.x,self.y
         elif self.show:
            plt.plot(self.x,self.y,'k-')
            plt.show()  
         else:
            pass
      else:
         print('>>> Analitical Solution is not available for this problem! <<<')

      #plt.plot(t,u,'-')
      #plt.show()
   
   def evalError(self, u : np.array , method : str , prints = False, write = False):
       self.x = np.linspace(self.t0,self.tf,self.n+1)
       self.y = np.array([self.u0 for i in range(self.n+1)])
       self.y = self.exactsolution(self.x,self.y).T

       if self.exactsolution == None: 
          print('Impossible evaluate the error, Analitical solutions is not provided')  
       #print(len(u), len(self.y))
       if len(u) != len(self.y):
          raise ValueError('Exception in EvalError : vector dimension must be equal')
       else:
               err = np.array([self.y[0] * i for i in range(len(self.y))])
       len(err)
       for i in range(len(self.y)):
            err[i] = np.abs(u[i] - self.y[i])
            
            if prints:
               if i==0 :
                    print('#----------- Errors of method %s ------- #' %method )
                    print()
               print('Error at point  %d :   %f' %(i,err[i]) )
       if write:
          filename = method + '.err'  
          with open(filename,'w') as f:
              for i in range(len(self.y)):
                   f.write('%.4f    %.6f\n' % (self.x[i],err[i]))
       
       maxErr = np.max(err) 
       print('The maximum error of the method %s  is :   %.6f' %(method,maxErr) )
       
 
       print('The mean error of the method %s  is :   %.6f' %(method, np.mean(err)) )
       return maxErr, err 
