#!/usr/bin/env python3

'''
   Module containing the Abstract Base Class for the hierarchy of the method to 
   solve System Of Differential Equation (SODE)

   
   The pourpose of this module is to make an interface for all the derived class 
   that directly implement (or group a branches of derived class) the solution 
   of SODE problem 
   
   @author Marco Ghiani, Glasgow 2018
'''

from abc import ABCMeta,abstractmethod

#--------------------------------------------------------------------------------------------

class Sode(metaclass=ABCMeta):
   '''  
        Interface class (Abstract Base) for all the sub classes that implement a method 
        to solve a System of Differential Equation
   '''
   @abstractmethod
   def solve(self,*args):    
        pass
   
   @abstractmethod
   def solve(self):    
      pass 
   
   def write2file(self):
       with open(self.file,'w') as f:
             for i in range(len(self.time)):
                f.write('%f  ' %self.time[i])
                for j in range(len(self.u[0])): 
                    f.write(' %f  ' %self.u[i,j])                    #, self.u[i,1], self.u[i,2]))
                f.write('\n')  

#--------------------------------------------------------------------------------------------






